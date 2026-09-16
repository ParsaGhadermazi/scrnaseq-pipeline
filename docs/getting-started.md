# Getting started

A first run, end to end, on data you can generate locally. Takes about ten
minutes once the image is pulled.

---

## 1. Requirements

| | |
|---|---|
| **Nextflow** | ≥ 24.04 — `curl -s https://get.nextflow.io \| bash` |
| **Docker** | running, with ~10 GB free |
| **Java** | 17+ (Nextflow needs it) |

Nothing else. Python, R, and every bioinformatics tool live inside the container.

```bash
docker pull parsaghadermazi/myrnaseqpipeline:0.2.0
```

---

## 2. Run the test suites first

This verifies your setup before you point the pipeline at real data, and it is
the fastest way to see what the pipeline actually does.

```bash
git clone https://github.com/ParsaGhadermazi/scrnaseq-pipeline.git
cd scrnaseq-pipeline

tests/run_tests.sh        /tmp/t1   # stage counts  — expect 11/11
tests/run_guard_tests.sh  /tmp/t2   # refusals      — expect 2/2
tests/run_qc_tests.sh     /tmp/t3   # stage qc      — expect 11/11
```

Each builds a synthetic dataset with **known ground truth** — a fixed number of
cells, genes, doublets and a planted batch effect — so results are checked by
equality rather than by eye. See [testing.md](testing.md).

---

## 3. The samplesheet

One file describes your experiment, and **every stage reuses it**. You never
re-enter the design, which is where samples get mismatched.

```csv
sample,srr,condition,batch
PBMC_healthy_1,SRR14575500,healthy,run1
PBMC_healthy_2,SRR14575501,healthy,run1
PBMC_treated_1,SRR14575502,treated,run2
PBMC_treated_2,SRR14575503,treated,run2
```

Only `sample` is required, plus a source (`srr`, or `fastq_1`/`fastq_2`).

| Column | Meaning |
|---|---|
| `sample` | unique id — output files are named from it |
| `srr` | SRA accession, downloaded automatically |
| `fastq_1`,`fastq_2` | local files, instead of `srr` |
| `condition` | **biological** group. Used by Part 5 (DE) |
| `batch` | **technical** group — chip, run, prep day. Defaults to `sample` |

Any other column you add rides along into `adata.obs` untouched — `donor`,
`timepoint`, `sex` all work with no code change.

> **`condition` and `batch` are not the same thing.** If every batch contains
> exactly one condition, they are the same variable and the pipeline will
> **refuse to run** — correcting the batch would delete your treatment effect.
> That is an experimental design problem, not a settings problem.

---

## 4. Stage 1 — reads to count matrices

```bash
nextflow run main.nf -profile docker \
  --stage counts --input samplesheet.csv --outdir results \
  --quantifier kallisto --genome_fasta genome.fa --gtf genes.gtf
```

If you do not pass `--index`, one is **built** and cached in
`--index_cache` (default `./references`) so it is built once, not per run.

**On Cell Ranger:** it is deliberately not in the container — 10x's licence
forbids redistribution. Supply it yourself and force amd64:

```bash
nextflow run main.nf -profile docker,amd \
  --stage counts --input samplesheet.csv \
  --quantifier cellranger \
  --cellranger_path /opt/cellranger-8.0.1 \
  --cellranger_ref  /refs/refdata-gex-GRCh38-2024-A
```

This needs ~64 GB RAM and will not run on a laptop.

**Output:** `results/h5ad/<sample>.h5ad` and `<sample>.raw.h5ad`.

> The `.raw.h5ad` is not a backup. The empty droplets in it **are** the ambient
> RNA measurement that stage `qc` needs. Keep it.

---

## 5. Stage 2 — QC through integration

```bash
nextflow run main.nf -profile docker \
  --stage qc --input samplesheet.csv --outdir results
```

It finds stage 1's output automatically. This runs empty-droplet detection,
ambient correction, doublet removal, cell and gene QC, then merges all samples
and corrects the batch effect.

**Output:** `results/qc/integrated/integrated.h5ad` — your analysis-ready object.

---

## 6. What you get

```
results/
├── h5ad/                       per-sample matrices (stage 1)
├── fastqc/  multiqc/           read QC
├── qc/
│   ├── emptydrops/             per-barcode calls
│   ├── ambient/                contamination rho per sample
│   ├── doublets/               per-cell doublet calls
│   ├── cell_qc/                what was filtered and why
│   └── integrated/
│       ├── integrated.h5ad     ← the finalized object
│       └── integration_metrics.json
└── pipeline_info/
    ├── run_manifest.txt        versions, params, container, git commit
    ├── software_versions.yml   captured from the running binaries
    └── timeline / report / trace / dag
```

Load it:

```python
import scanpy as sc
a = sc.read_h5ad("results/qc/integrated/integrated.h5ad")

a.layers["counts"]       # raw counts — use these for DE
a.obsm["X_integrated"]   # batch-corrected — use for clustering/UMAP
```

---

## Where to go next

| | |
|---|---|
| [architecture.md](architecture.md) | diagrams of both stages, and why they are shaped that way |
| [parameters.md](parameters.md) | every parameter and default |
| [usage.md](usage.md) | detailed per-stage usage |
| [lectures/](lectures/) | the biology and methods, worked through from first principles |
| [COVERAGE.md](COVERAGE.md) | every tutorial tool vs what is implemented, and what is verified |
