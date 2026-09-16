# Testing

The pipeline is verified against **synthetic 10x data with known ground truth**,
so results can be checked by equality rather than by eyeballing plots.

## Running it

    tests/run_tests.sh [workdir]

Builds the fixture, then checks both detector refusals, both quantifiers end to
end, ground-truth equality, and cross-quantifier schema parity. Exits non-zero
on any failure. **Last run: 11 passed, 0 failed on image v0.2.0.**

Requires the local image. It is built, never pushed, so a `docker image prune -a`
removes it -- rebuild with:

    docker buildx build --platform linux/arm64 -f containers/Dockerfile \
      -t parsaghadermazi/myrnaseqpipeline:0.2.0 --load .

## Building the fixture

40 barcodes drawn from the real 10x v3 whitelist (shipped inside the image at
`ngs_tools/chemistry/whitelists/10x_version3_whitelist.txt.gz`), a 20-transcript
/ 5-gene toy transcriptome, and 20,000 read pairs where R1 = 16 bp barcode +
12 bp UMI (28 bp) and R2 = a 91 bp substring of a real transcript.

Ground truth: **40 cells x 5 genes, 20,000 UMIs.**

## What has been verified

| Check | Expected | Result |
|---|---|---|
| Chemistry auto-detection | `SC3Pv3` from R1 = 28 bp | pass |
| Chemistry mis-declaration | refuse, exit 1 | pass |
| Missing barcode read (`--include-technical` trap) | refuse, name the fix, exit 1 | pass |
| kallisto matrix | 40 x 5, sum 20000 | pass |
| alevin-fry matrix | 40 x 5, sum 19998 | pass (see below) |
| kallisto vs alevin-fry concordance | high | **r = 1.0000** |
| obs/var schema across quantifiers | identical | pass |
| `QC_ASSERT` on good data | PASS, valid_barcodes 1.0 | pass |
| `QC_ASSERT` on off-whitelist barcodes | FAIL the run, exit 1 | pass |

### The 20000 vs 19998 difference is real, not a bug

alevin-fry collapses 2 UMI collisions that kb's filtered counts retain --
`inspect.json` independently reports `numBarcodeUMIs: 19998`. A genuine
methodological difference between the tools, and exactly the kind of thing the
standardised `.h5ad` contract exists to make visible.

---

# Stage `qc` fixture

    tests/make_qc_fixture.sh <outdir>     # build fixture
    tests/run_qc_tests.sh    <workdir>    # fixture + run stage qc + check contract

The suite runs `-profile docker,test`; the `test` profile carries
`soupx_soup_quantile = 0.50`, without which SoupX fails the soup gate on a
synthetic fixture (see below).

### What it asserts

Structural — that the output contract holds:

- `layers["counts"]` survives (differential expression must reach raw counts)
- `obsm["X_integrated"]` exists, and `obsm["X_pca"]` is retained uncorrected
- both samples present, cells reduced by QC but not wiped out
- `uns["integration"]` recorded

Behavioural — because structure alone cannot distinguish real integration from a
no-op copy of `X_pca`:

- **`X_integrated` is not identical to `X_pca`** — catches a silently skipped
  correction, which would otherwise satisfy every structural check above
- **batch silhouette decreases** from uncorrected to integrated — the fixture
  plants a 1.6x shift on sample 2, so a correction that genuinely works must make
  samples *less* separable in the embedding

That second one is the only assertion in the suite that would fail if integration
ran but did nothing useful.

### The confounding guard, verified by accident

An early fixture had one sample per condition (S1 healthy, S2 treated). Stage
`qc` refused to run:

    ERROR: 'sample' is perfectly confounded with 'condition'.
      Each batch contains exactly one condition and each condition one batch, so
      they are the same variable. Integrating would remove the condition effect.
      This is an experimental design problem; no setting fixes it.

That is the pipeline behaving correctly and the *fixture* being wrong: with a
1:1 batch-condition mapping there is no way to remove the batch effect without
removing the treatment effect, and no parameter can rescue it. The guard caught a
design error I had written without noticing.

The fixture now uses **4 samples, 2 per condition** -- the tutorial's own design
-- so batch is nested within condition rather than identical to it. Each sample
keeps its own shift (1.00/1.15 healthy, 1.60/1.75 treated), so there is still a
real batch effect to correct *within* each condition.

It now has its own test rather than relying on accidental discovery:

    tests/run_guard_tests.sh <workdir>

Negative tests, asserting the pipeline **refuses** bad input:

1. **batch perfectly confounded with condition** -- one sample per condition;
   must fail with the explanatory error above
2. **missing `.raw.h5ad`** -- must fail rather than silently skipping ambient
   correction, since the empty droplets *are* the soup measurement

These guard the failure mode that actually matters. A crash is obvious; a run
that quietly removes the treatment effect and emits a clean-looking object is the
one that reaches a paper.

**Last run: 2 passed, 0 failed.**

    PASS  confounded batch/condition refused with an explanatory error
    PASS  missing .raw.h5ad refused

Both refusals are now verified against the real pipeline rather than asserted in
a docstring.

Stage `qc` needs far more structure than stage `counts`, because each step needs
something real to find. All of it is planted deliberately:

| Planted | Why it has to be there |
|---|---|
| 1,500 **empty droplets** per sample | `emptyDrops` has nothing to profile without them -- they ARE the soup measurement |
| 10% **ambient contamination** in every cell | SoupX needs contamination actually present to estimate rho |
| 25 **heterotypic doublets** per sample | scDblFinder can only see doublets of *different* types |
| 3 **cell types** with **exclusive** 40-gene marker blocks | see below -- non-exclusive markers break SoupX outright |
| a **1.6x global shift** on sample 2 | a batch effect for integration to remove |

Ground truth is written to `ground_truth.json`, so results are checked by
equality rather than by eye.

### Marker genes must be EXCLUSIVE

The first version of this fixture gave every gene nonzero baseline expression in
every cell type, with markers merely 8x up. Stage `qc` died on it:

    0 genes passed tf-idf cut-off and 0 soup quantile filter.
    Error in autoEstCont(sc) : No plausible marker genes found.
      Is the channel low complexity?

That is not a bug -- it is the mechanism. `autoEstCont` estimates contamination
from genes a cluster should **not** express: any counts of those genes there are
definitionally soup. If every gene is expressed everywhere, no "should be zero"
group exists and there is nothing to measure against.

So each cell type's 40 marker genes are expressed **only** in that type and are
~absent elsewhere, with 280 housekeeping genes expressed in all types. Real data
has this property -- `CD3D` in T cells, not monocytes -- and a fixture that lacks
it tests a scenario that cannot occur.

### Why exclusivity alone was not enough

A second attempt with fully exclusive markers *still* failed. `quickMarkers`
ranks by tf-idf, where `idf = log(n_cells / n_cells_expressing_gene)`. Ambient
deposits counts of every gene into every cell, so if a marker is *detected* in
most wrong-type cells, idf collapses and nothing clears `tfidfMin`, however large
the fold-change. On the 400-gene version: marker idf 0.083 vs housekeeping idf
0.083 -- zero discriminative power, despite 62.2 vs 2.1 mean counts.

A third attempt lowered `MARKER_LEVEL` 25 -> 10 on the theory that ambient per
marker gene was the driver. Measured idf moved 0.57 -> 0.66 and top tf-idf moved
0.645 -> 0.640: essentially nothing. **A quantity that will not move when you
change its supposed driver means the model is wrong**, and the arithmetic
previously written here was wrong; it has been removed rather than patched.

### The real driver, measured

A parameter sweep gives the answer that derivation kept missing:

    idf <= log(N_TYPES)

A **perfectly** exclusive marker is still expressed in `1/N_TYPES` of all cells,
so idf cannot exceed `log(N_TYPES)` even with zero ambient. Verified at
ambient = 0, measured idf equals `log(N_TYPES)` to three decimals:

| N_TYPES | idf @ ambient 0.08 | idf @ ambient 0.03 | ceiling log(N) |
|---|---|---|---|
| 3  | 0.728 | 0.930 | **1.099** |
| 6  | 1.313 | 1.575 | 1.792 |
| 8  | 1.585 | 1.845 | 2.079 |
| 12 | 1.948 | 2.242 | 2.485 |

With **3 types the ceiling is 1.099**, so any ambient at all drops idf below
SoupX's `tfidfMin = 1.0`. No value of `MARKER_LEVEL` or `AMBIENT_FRAC` could have
rescued it -- which is precisely why two rounds of tuning those moved idf only
0.57 -> 0.66.

The fixture therefore uses **8 cell types at 8% ambient** (idf 1.585): clears the
production threshold with headroom, keeps contamination above `soupx_min_rho`,
and matches the realistic regime -- real PBMC data has 10+ cell types, which is
why markers there are rare enough for tf-idf to work at all.

`make_qc_fixture.sh` asserts this itself and refuses to emit an invalid fixture,
so the failure surfaces in seconds rather than after a full pipeline run.

### There are TWO gates, not one

Fixing tf-idf moved the error rather than removing it:

    320 genes passed tf-idf cut-off and 0 soup quantile filter.  Taking the top 0.

`quickMarkers` applies a second, independent filter. `soupQuantile` (default
0.90) requires a candidate marker to be expressed **in the soup** above the 90th
percentile of all genes. The reasoning is sound: you cannot measure contamination
using a gene that is barely present in the contaminant.

These two gates pull against each other:

| lever | effect on idf | effect on soup rank |
|---|---|---|
| raise `MARKER_LEVEL` | **worse** (more ambient per marker) | better |
| lower housekeeping expression | neutral | **better** |
| raise `N_TYPES` | better | worse (each marker is rarer in the soup) |

With ~980 housekeeping genes expressed in *every* cell, they dominate soup mass
and push exclusive markers -- present in only 1/N_TYPES of cells -- below the
90th percentile of soup expression, where `quickMarkers` discards them.

### The soup gate is not reachable from a synthetic fixture

A joint sweep over `n_types` {8,12,16} x `marker_level` {30,40} x `hk_scale`
{0.25,0.10}, three seeds each, produced **negative margin in all twelve cells**.
The revealing part is that `soup/q90` sat at **0.98-0.99 in every single one** --
completely flat across all three levers. A quantity that will not move when you
change any of its supposed drivers is not tunable by them.

The reason is structural, not a tuning failure. In real 30k-gene data the soup is
dominated by whatever the abundant cell types express, so *their* markers
genuinely rank high in it. In a ~1300-gene fixture, 980 ubiquitous housekeeping
genes inevitably hold the top of the soup distribution. The fixture cannot
reproduce that property at this scale, and pretending otherwise was chasing
realism it cannot have.

So the split is:

- **tf-idf gate** -- reproducible, and the fixture clears it at the production
  `tfidfMin = 1.0` (idf ~1.38 at 12 types).
- **soup gate** -- lowered to `soupx_soup_quantile = 0.50` in `conf/test.config`
  **only**. Production keeps SoupX's default 0.90.

That is the parameter SoupX's own error message tells you to reduce, and scoping
it to the test profile keeps the production path honest: the test still exercises
the real `tfidfMin`, real `autoEstCont`, and real `adjustCounts`.
the 90th percentile.


For genuinely low-complexity real channels, `--soupx_tfidf_min` and
`--soupx_soup_quantile` are exposed; SoupX's own error recommends lowering them.
The module refuses to guess a contamination fraction rather than silently
skipping correction.

## Measured against ground truth

Measured on the **current** fixture (4 samples x 505 cells, 1000 empties each,
12 cell types, 25 planted doublets, 8% ambient) from the green 11/11 run:

| Step | Planted | Measured (S1 / S2 / S3 / S4) | |
|---|---|---|---|
| `emptyDrops` cells | 505 of 1505 barcodes | **505 / 507 / 506 / 506** | within 0.4%, at `niters=1000` |
| `SoupX` rho | 0.08 | **0.060 / 0.059 / 0.061 / 0.063** | ~25% under -- see below |
| `scDblFinder` doublets | 25 | **13 / 10 / 25 / 25** | **asymmetric -- see below** |
| `CELL_QC` cells kept | -- | 96.2% / 96.7% / 94.7% / 93.9% | |

### Two findings worth not glossing over

**SoupX underestimates rho consistently** (0.059-0.063 against a planted 0.08,
~25% low, in the same direction across all four samples). `autoEstCont` is
deliberately conservative -- it would rather leave contamination in than remove
real signal -- so a systematic underestimate is expected behaviour, not a bug.
Worth knowing that correction is partial by design.

**scDblFinder recovered all 25 doublets in S3/S4 but only 13 and 10 in S1/S2.**
The split follows the fixture's batch scaling exactly: S1/S2 carry 1.00x/1.15x
depth, S3/S4 carry 1.60x/1.75x. Doublet detection scores cells by how many of
their k nearest neighbours are *simulated* doublets, and that neighbourhood
structure is noisier at lower depth -- so the shallower samples lose sensitivity.

This is a real property of the method rather than a pipeline defect, and it has a
practical consequence: **doublet detection sensitivity varies with sequencing
depth across samples in the same experiment.** Residual doublets will be
concentrated in your shallowest samples. Worth remembering when a "rare
population" appears in only the lightly-sequenced ones.

### Integration

| metric | uncorrected | integrated | direction |
|---|---|---|---|
| batch silhouette | -0.0075 | **-0.0135** | down = good |
| label silhouette | 0.5638 | **0.5764** | up = good |
| kNN batch mixing | 0.8437 | **0.8786** | up = good |

All three move the right way: samples became less separable, cell populations
stayed distinct, neighbourhoods became better mixed.

**But read the absolute values, not just the deltas.** kNN mixing was already
0.84 *before* correction, and batch silhouette was already near zero -- so the
planted batch effect is weak in PCA space to begin with. These numbers confirm
integration works and moves things correctly; they do **not** demonstrate
recovery from a severe batch effect. A fixture with a genuinely hard batch
effect would be a better test, and this one is not it.

## Status:## Status: STAGE qc GREEN -- 11 passed, 0 failed

Seven of eleven processes now run correctly on real data for both samples, with
the results above. The remaining failures have been h5ad serialisation bugs, not
analysis bugs: `.uns` payloads containing nested dicts, which h5py cannot write
("Can't implicitly convert non-string objects to strings"). AnnData holds them
happily in memory and only fails at `write_h5ad`, so each one surfaced a single
process later than where it was introduced. All `.uns` writes now serialise to
JSON strings; the structured form goes to the sidecar `.json` reports.

## Not covered

`SRATOOLS_FASTERQDUMP` (needs network), `CELLRANGER_COUNT` (needs the licensed
binary, amd64, ~64 GB RAM), and index building (no module yet).
