# Parameters

Every parameter and its default, generated from `nextflow.config`. Override any
of them on the command line (`--param value`) or in a config file.

> Defaults are chosen to **generalise**, not to reproduce the tutorial exactly.
> Where the two differ it is called out below. `--qc_method fixed` reproduces the
> tutorial's thresholds.

---

## Core

| Parameter | Default | Notes |
|---|---|---|
| `--input` | *required* | samplesheet CSV — see [schema](../assets/schema_input.json) |
| `--outdir` | `results` | |
| `--stage` | `counts` | `counts` \| `qc` \| `annotate` \| `de` |
| `--step` | `fastq` | within `counts`: `fastq` starts from reads, `counts` from matrices |
| `--container_image` | `parsaghadermazi/myrnaseqpipeline:0.2.0` | changed in one place only |
| `--max_cpus` / `--max_memory` / `--max_time` | `8` / `14.GB` / `48.h` | per-process ceilings |

## Stage `counts` — quantification

| Parameter | Default | Notes |
|---|---|---|
| `--quantifier` | `cellranger` | `cellranger` \| `alevinfry` \| `kallisto`. One per run |
| `--cellranger_path` | `null` | **required for cellranger** — bind-mounted, never containerised (10x licence) |
| `--cellranger_ref` | `null` | prebuilt `refdata-gex-*` dir |
| `--index` | `null` | prebuilt index; if absent one is **built** from genome + GTF |
| `--genome_fasta` / `--gtf` | `null` | needed only when building an index |
| `--chemistry` | `null` | **auto-detected from R1 length.** A declaration contradicting the read length is a hard error |
| `--expect_cells` | `null` | soft prior only since Cell Ranger v3; `--force-cells` is the dangerous one |
| `--skip_fastqc` / `--skip_multiqc` | `false` | |

## Reference / index building

| Parameter | Default | Notes |
|---|---|---|
| `--index_cache` | `<launchDir>/references` | `storeDir` cache — built once, reused. **Point at an external drive if disk is tight** |
| `--reference_id` | `null` | derived from fasta + GTF names |
| `--gtf_biotypes` | protein_coding, lncRNA, IG/TR genes | filtering is **not cosmetic** — a pseudogene overlapping its parent silently deletes counts from the real gene |
| `--af_ref_type` | `spliced+intronic` | **splici**. Keeps spliced/unspliced separable, which is what makes RNA velocity (Part 13) possible at all |
| `--read_length` | `91` | R2 length; sets splici intronic flank |
| `--kb_workflow` | `standard` | `nac` is kallisto's velocity-capable equivalent |

## Stage `qc` — empty droplets

| Parameter | Default | Notes |
|---|---|---|
| `--counts_dir` | `<outdir>/h5ad` | where stage `counts` left its output |
| `--emptydrops_skip` | `false` | Cell Ranger already applies a similar call |
| `--emptydrops_lower` | `100` | barcodes below this **define** the ambient profile |
| `--emptydrops_niters` | `10000` | Monte Carlo; the test profile drops this to 1000 for speed |
| `--emptydrops_fdr` | `0.01` | |
| `--seed` | `42` | emptyDrops and scDblFinder are stochastic |

## Stage `qc` — ambient RNA

| Parameter | Default | Notes |
|---|---|---|
| `--ambient_method` | `soupx` | `soupx` \| `decontx`* \| `cellbender`* \| `none`. **Single value** — it rewrites the whole matrix |
| `--soupx_min_rho` | `0.05` | correct only if contamination ≥ 5% (tutorial's rule) |
| `--soupx_rho` | `null` | bypass `autoEstCont` with a known fraction |
| `--soupx_tfidf_min` | `1.0` | lower if `autoEstCont` finds no marker genes |
| `--soupx_soup_quantile` | `0.90` | lower for low-complexity channels |

<sub>*planned, not yet implemented — see [COVERAGE.md](COVERAGE.md)</sub>

## Stage `qc` — doublets

| Parameter | Default | Notes |
|---|---|---|
| `--doublet_method` | `scdblfinder` | **a list is allowed** — each method emits only a CSV, so running several is nearly free |
| `--doublet_removal` | `primary` | `primary` \| `consensus` (needs ≥2 methods) \| `none` |
| `--doublet_rate` | `null` | estimated from cell count (~0.8% per 1,000) |

## Stage `qc` — cell and gene filtering

| Parameter | Default | Notes |
|---|---|---|
| `--qc_method` | `mad` | `mad` generalises across tissue and depth; `fixed` reproduces the tutorial |
| `--mad_nmads` | `3` | |
| `--max_mito_pct` | `20` | **absolute ceiling, always applied** alongside MAD |
| `--min_cells_per_gene` | `3` | |
| `--remove_hb_genes` | `true` | **set `false`** for bone marrow / erythropoiesis, where haemoglobin genes *are* the biology |
| `--min_genes` / `--max_genes` | `500` / `5000` | `fixed` mode only |
| `--min_counts` / `--max_counts` | `800` / `20000` | `fixed` mode only |

## Stage `qc` — normalisation and integration

| Parameter | Default | Notes |
|---|---|---|
| `--norm_method` | `lognorm` | `lognorm` \| `pearson_residuals` |
| `--n_hvg` | `2000` | selected **batch-aware**, or donor and stress genes dominate |
| `--n_pcs` / `--n_pcs_merged` | `30` / `50` | |
| `--prelim_resolution` | `0.8` | throwaway clustering for SoupX only |
| `--integration_method` | `harmony` | `harmony` \| `scanorama` \| `none`. `cca`/`rpca`/`fastmnn` are planned R modules |
| `--batch_key` | `sample` | **technical** grouping — never `condition` |
| `--condition_key` | `condition` | used only by the confounding guard |
| `--mixing_k` | `50` | neighbours for the kNN mixing metric |
| `--silhouette_subsample` | `5000` | silhouette is O(n²) |

## Correctness gates

| Parameter | Default | Notes |
|---|---|---|
| `--qc_strict` | `true` | fail the run rather than warn |
| `--min_valid_barcode_fraction` | `0.70` | below this signals the wrong chemistry |
| `--fractional_counts_policy` | `round` | alevin-fry's EM emits fractional counts; DESeq2 and scDblFinder assume integers |

---

## Profiles

```bash
-profile docker          # pull the published image
-profile docker,test     # + small/fast settings for the test fixtures
-profile conda           # build environments from env/*.yml instead
-profile amd             # force linux/amd64 (required for Cell Ranger)
-profile arm             # force linux/arm64
-profile slurm           # SLURM executor
```

`test` lowers `emptydrops_niters` to 1000, `silhouette_subsample` to 2000, and
`soupx_soup_quantile` to 0.50 — the last because SoupX's soup-quantile gate is
structurally unreachable on a small synthetic fixture. Production keeps 0.90.
See [testing.md](testing.md).
