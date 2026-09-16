# Usage

The pipeline runs **one stage at a time**. Each stage is a self-contained
invocation that finds the previous stage's output by convention.

| Stage | Tutorial parts | Input | Output |
|---|---|---|---|
| `counts` | 1 | samplesheet + FASTQ | per-sample `.h5ad` + `.raw.h5ad` |
| `qc` | 2, 2-2, 3 | per-sample `.h5ad` | **one merged, integrated `.h5ad`** |
| `annotate` | 4 | merged `.h5ad` | cell-type labels *(not built)* |
| `de` | 5 | labelled `.h5ad` | DE / pathways *(not built)* |

## The samplesheet is the experiment manifest

Declared **once** and reused by every stage, so experimental design is never
re-entered (which is where samples get mismatched).

```csv
sample,srr,condition,batch
PBMC_healthy_1,SRR14575500,healthy,run1
PBMC_treated_1,SRR14575501,treated,run1
```

Only `sample` is required, plus either `srr` or `fastq_1`/`fastq_2`.
Any column that isn't reserved is carried verbatim into the meta map and lands
in `adata.obs` — so `donor`, `timepoint`, `sex` work with no code change.

`batch` is the **technical** grouping (chip, run, prep day), defaults to
`sample`, and is consumed by integration. It is never the same thing as
`condition`.

## Stage 1 — counts

```bash
nextflow run main.nf -profile docker \
  --stage counts --input samplesheet.csv --outdir results \
  --quantifier kallisto --genome_fasta genome.fa --gtf genes.gtf
```

`--quantifier` is one of `cellranger | alevinfry | kallisto`. If you don't pass
an index, one is **built** from `--genome_fasta` + `--gtf` and cached in
`--index_cache` (default `<launchDir>/references`) so it is built once, not
per run. See [reference.md](reference.md).

Cell Ranger is never containerised — 10x's licence forbids redistribution — so
it needs `--cellranger_path` pointing at your unpacked tarball, and `-profile amd`.

## Stage 2 — qc

```bash
nextflow run main.nf -profile docker \
  --stage qc --input samplesheet.csv --outdir results
```

Finds `results/h5ad/<sample>.h5ad` and `<sample>.raw.h5ad` automatically;
override with `--counts_dir`.

> **The `.raw.h5ad` is not optional.** `emptyDrops` and SoupX measure the ambient
> soup *from the empty droplets*, which exist only in the raw matrix. The stage
> refuses to start without it.

Runs: `emptyDrops` → preliminary clustering → SoupX → scDblFinder → cell/gene QC
→ merge → normalise/HVG → integrate → integration QC.

### Parameters worth knowing

| Param | Default | Note |
|---|---|---|
| `--qc_method` | `mad` | `mad` generalises; `fixed` reproduces the tutorial exactly |
| `--mad_nmads` | `3` | |
| `--max_mito_pct` | `20` | absolute ceiling, **always** applied alongside MAD |
| `--ambient_method` | `soupx` | single value — it rewrites the whole matrix |
| `--soupx_min_rho` | `0.05` | skip correction below this contamination fraction |
| `--soupx_tfidf_min` | `1.0` | lower if `autoEstCont` reports "no plausible marker genes" |
| `--soupx_soup_quantile` | `0.90` | lower for low-complexity channels; see note below |
| `--soupx_rho` | `null` | bypass `autoEstCont` with a known contamination fraction |
| `--doublet_method` | `scdblfinder` | a list is allowed; each emits only a CSV |
| `--doublet_removal` | `primary` | `consensus` needs ≥2 methods |
| `--remove_hb_genes` | `true` | **set `false`** for bone marrow / erythropoiesis, where haemoglobin genes *are* the biology |
| `--integration_method` | `harmony` | corrects an embedding, never expression |
| `--batch_key` | `sample` | technical grouping, never `condition` |

**If SoupX fails with "No plausible marker genes found"**, that is not a crash to
work around — it means `autoEstCont` could not find genes a cluster *should not*
express, which is the only thing contamination can be measured against. Either
the channel genuinely is low-complexity (lower `--soupx_tfidf_min` /
`--soupx_soup_quantile`), or supply `--soupx_rho` from a known estimate. The
module refuses to guess a contamination fraction rather than silently skipping
correction and emitting uncorrected counts labelled as corrected.

Two guards will stop the run rather than produce quietly wrong output:

- **one batch level** → integration is skipped and says so
- **`batch_key` perfectly confounded with `condition_key`** → refuses, because
  correcting the batch would delete the treatment effect. That is an
  experimental-design problem no setting can fix.

## Output contract

```
layers["counts"]      raw counts, untouched      -> differential expression
X                     log-normalised             -> marker plots
obsm["X_pca"]         uncorrected
obsm["X_integrated"]  batch-corrected            -> neighbours, clustering, UMAP
```

Nothing downstream ever reads "corrected expression" — it does not exist.
Integration modifies the embedding only.

## Provenance

Every run writes `results/pipeline_info/`: `run_manifest.txt` (versions, params,
container, git commit, timings), `software_versions.yml` (versions captured from
the running binaries), plus timeline, report, trace and DAG.

## Verification status

**Stage `counts` is tested** — 11/11 checks, see [testing.md](testing.md).

**Stage `qc` is tested** — 11/11 checks on image v0.2.0, all 32 tasks succeeding.
The refusal guards pass 2/2 (`tests/run_guard_tests.sh`).

**Cell Ranger and the three index builders have genuinely never executed** --
they need the licensed x86_64 binary and ~32-64 GB RAM, so they are deferred to
cluster testing.

See [COVERAGE.md](COVERAGE.md) for per-module status and [testing.md](testing.md)
for what each suite actually asserts.
