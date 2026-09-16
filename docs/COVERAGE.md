# Tutorial coverage

Tracks every tool and method in the ngs101 single-cell series against this pipeline.

**Standing rule (revised 2026-09-09):** every tool the tutorial uses gets
implemented. A cross-language alternative is added **only when it is genuinely
competitive**, not for symmetry -- if the best tool for a step exists in one
language only, we use that language and say so. Where no real counterpart
exists, this file states it rather than skipping quietly.

Consequence: **Part 2 is R-led**, because the best ambient-RNA and doublet tools
are R-only. That is a deliberate choice, not an oversight.

Status vocabulary:

- `verified` -- executed on real data and checked against planted ground truth
- `runs` -- executes without error; no ground-truth assertion yet
- `written, unrun` -- code exists, has never executed
- `planned` -- not written
- `n/a` -- no counterpart exists

A figure in this table is only meaningful alongside the fixture **and the run**
that produced it. Figures from ad-hoc probes are not recorded here -- only ones a
test suite produced. When the fixture changes, figures are retired rather than
carried forward.

**Current execution status** (image `parsaghadermazi/myrnaseqpipeline:0.2.0`, arm64):

- stage `counts` -- **11/11 passing**, cross-quantifier concordance r = 1.0000
- refusal guards -- **2/2 passing** (`tests/run_guard_tests.sh`)
- stage `qc` -- **11/11 passing**, 32 tasks, 8m32s
- Cell Ranger modules and the three index builders -- **never executed**
  (need the licensed x86_64 binary and ~32-64 GB RAM; deferred to cluster testing)

---

## Part 1 -- FASTQ to count matrix

| Step | Tutorial tool | Python / alt counterpart | Module | Status |
|---|---|---|---|---|
| Download from SRA | `prefetch` + `fasterq-dump` | none needed | `SRATOOLS_FASTERQDUMP` | written, unrun |
| Identify read roles | manual inspection | **`detect_10x_reads.py`** (ours) | `TENX_READ_DETECT` | **tested** |
| Rename to 10x convention | manual `mv` | folded into above | `TENX_READ_DETECT` | **tested** |
| Read QC (R2 only) | `fastqc` | none | `FASTQC` | **tested** |
| Aggregate QC | `multiqc` | (is Python) | `MULTIQC` | **tested** |
| Build reference/index | `cellranger mkgtf` + `mkref` | `simpleaf index` (splici), `kb ref` | `CELLRANGER_MKREF`, `SIMPLEAF_INDEX`, `KB_REF` | written, **unrun** |
| Quantification | `cellranger count` | `alevin-fry`/`simpleaf`, `kallisto\|bustools` | `CELLRANGER_COUNT` | written, unrun |
| " | " | " | `SIMPLEAF_QUANT` | **tested** |
| " | " | " | `KB_COUNT` | **tested** |
| Run metric assertions | eyeball `web_summary.html` | **`QC_ASSERT`** (ours) | `QC_ASSERT` | **tested** (pass + fail) |
| Standardise output | n/a | `anndata` | `COUNTS_TO_H5AD` | **tested** (both quantifiers) |

Verified by `tests/run_tests.sh`: **11 passed, 0 failed on image v0.2.0**
(cross-quantifier concordance r = 1.0000).

### Known gaps in Part 1

**Decision (2026-09-01): these are deferred to cluster testing, not blockers.**
Part 1 is considered complete. The modules below are written against documented
interfaces and will be debugged on first real use on an HPC system.


- **Index building is implemented but never executed.** If no index is supplied
  the pipeline builds one from `--genome_fasta` + `--gtf`, cached via `storeDir`.
  All three modules are written against documented interfaces and unrun --
  `CELLRANGER_MKREF` needs the licensed binary and ~32 GB+ RAM, and the splici
  path of `SIMPLEAF_INDEX` has not been exercised. See `docs/reference.md`.
- `SRATOOLS_FASTERQDUMP` has never run: needs network and a real accession.
- `CELLRANGER_COUNT` has never run: needs the licensed binary, amd64, and ~64 GB RAM.

Deliberately excluded: **STARsolo** -- same family as Cell Ranger, near-identical
output, no additional failure mode covered.

Not in the tutorial, added by us: `TENX_READ_DETECT` (kills the `--include-technical`
trap, SRA file-ordering ambiguity, and silent chemistry mis-declaration) and
`QC_ASSERT` (fails loudly on a low valid-barcode fraction instead of passing a
garbage matrix downstream).

## Part 2 -- QC and cell filtering

> **Architecture constraint:** R modules can READ `.h5ad` but cannot WRITE it.
> Every R step here emits CSV annotations; only ambient correction returns a
> matrix, and it does so as MTX for a Python step to convert. See
> [python-r-bridge.md](python-r-bridge.md) -- this shapes every module below.

| Step | Tutorial tool | Python counterpart | Status |
|---|---|---|---|
| Empty droplet detection | `DropletUtils::emptyDrops` | CellBender | `EMPTYDROPS` | **verified** on an earlier fixture; re-measure pending |
| Ambient RNA correction | `SoupX` | CellBender | `AMBIENT_SOUPX` | **verified** on an earlier fixture; re-measure pending |
| " (alternative) | `DecontX` (celda) | -- | -- | planned |
| Doublet detection | `scDblFinder` | Scrublet, DoubletDetection | `DOUBLETS_SCDBLFINDER` | **verified** on an earlier fixture; re-measure pending |
| " (alternative) | `DoubletFinder` | -- | -- | planned |
| Cell/gene QC + filtering | Seurat `subset()`, fixed thresholds | scanpy | `CELL_QC` | **verified** |
| QC thresholds -- adaptive | `scuttle::isOutlier` (MAD) | own impl, MAD==0 guarded | `CELL_QC --method mad` | written, unrun |
| QC thresholds -- model-based | `miQC` | none | -- | planned |
| Normalisation + HVG | Seurat | scanpy, batch-aware | `NORMALIZE_HVG` | **verified** |

## Part 3 -- integration

| Step | Tutorial tool | Python counterpart | Module | Status |
|---|---|---|---|---|
| Merge samples | Seurat split layers | `anndata.concat` (inner join) | `MERGE_SAMPLES` | **verified** |
| Integration | Harmony | `harmonypy` called directly (see note) | `INTEGRATE` | **verified** by the qc suite (batch silhouette -0.0075 -> -0.0135) |
| " (alternatives) | CCA, RPCA, FastMNN | scanorama | -- | planned (R module) |
| Confounding guard | -- | ours: refuses batch == condition | `INTEGRATE` | **verified** (`run_guard_tests.sh` 2/2) |
| Silhouette scoring | `cluster::silhouette` | `sklearn` | `INTEGRATION_QC` | runs; emits metrics |
| Local mixing | LISI | ours: kNN entropy | `INTEGRATION_QC` | runs; emits metrics |

**Note on Harmony:** `sc.external.pp.harmony_integrate` is stale relative to
harmonypy 2.0.0 -- it transposes an already-correctly-oriented matrix. We call
`harmonypy.run_harmony` directly and orient on `n_obs`. See
[python-r-bridge.md](python-r-bridge.md#upstream-version-traps).

Clustering is deliberately NOT here -- it opens stage `annotate`, so re-clustering
at a new resolution does not re-run ambient correction, doublets and integration.

**Note on thresholds.** The tutorial's main workflow uses manual, visually-set
fixed thresholds (nFeature 500-5000, nCount 800-20000, percent.mt <10%) and
argues for visual inspection over "automatic statistical methods". It gives
tissue-specific guidance (PBMC <10%, brain <5%, tumor <20%) and, in a reader
reply about scaling to 20-40 samples, points to `scuttle::isOutlier` (MAD) and
`miQC` as the right approach at scale.

We implement both: `qc_method = 'fixed'` reproduces the tutorial exactly,
`qc_method = 'mad'` is the default for multi-sample runs.

## Part 2-2 -- QC tool comparison
SoupX vs DecontX; DoubletFinder vs scDblFinder. All four required. | planned

## Parts 3-17
Integration/clustering, annotation, DE + pathways, object internals, Monocle 3,
Slingshot, PDX, CellChat, NicheNet, CopyKAT, hdWGCNA, scVelo, CellRank,
scplotter, decoupleR x2. Detail filled in as each lecture is reached.

**Known R-only (no Python counterpart worth the name):** CellChat, NicheNet, hdWGCNA.
**Known Python-only (no R counterpart):** scVelo, CellTypist.

## scATAC-seq track
Parts 1-3 (Signac). In scope, scheduled after the scRNA-seq track completes.
