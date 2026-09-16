# myRNASeqPipeline

A single-cell RNA-seq pipeline built step by step against the
[ngs101 tutorial series](https://ngs101.com/tutorials/#single-cell-seq), with the
reasoning written down as it was worked out.

Nextflow DSL2 orchestration over Python and R scripts, in one multi-arch
container. Runs **one stage at a time**.

## Status

| Stage | Parts | State |
|---|---|---|
| `counts` — FASTQ → count matrices | 1 | **tested**, 11/11 on image v0.2.0 |
| `qc` — QC, ambient, doublets, merge, integrate | 2, 2-2, 3 | **tested**, 11/11 on image v0.2.0 |
| `annotate` | 4 | not started |
| `de` | 5 | not started |

Both stages run end to end against synthetic fixtures with known ground truth,
and the refusal guards pass. Stage `qc`'s suite asserts not just that integration
*ran* but that it **did something**: `X_integrated` must differ from `X_pca`, and
batch silhouette must decrease against the planted batch effect.

Known unverified: Cell Ranger modules (need the licensed x86_64 binary and
~32–64 GB RAM) and the three index builders — deferred to cluster testing. See
[`docs/COVERAGE.md`](docs/COVERAGE.md) for per-module status; nothing is marked
done that has not run, and figures are retired when the fixture they were
measured on changes.

## Quickstart

```bash
# stage 1 — reads to count matrices
nextflow run main.nf -profile docker \
  --stage counts --input samplesheet.csv --outdir results \
  --quantifier kallisto --genome_fasta genome.fa --gtf genes.gtf

# stage 2 — QC through integration
nextflow run main.nf -profile docker \
  --stage qc --input samplesheet.csv --outdir results
```

Samplesheet — only `sample` plus a source is required:

```csv
sample,srr,condition,batch
PBMC_healthy_1,SRR14575500,healthy,run1
```

## Container

Published multi-arch on Docker Hub — `linux/amd64` and `linux/arm64`:

```bash
docker pull parsaghadermazi/myrnaseqpipeline:0.2.0
```

`params.container_image` already points at it, so `-profile docker` pulls it
automatically. To rebuild locally instead:

```bash
docker buildx create --name scrnabuilder --driver docker-container --bootstrap   # once
docker buildx build --builder scrnabuilder --platform linux/amd64,linux/arm64 \
  -f containers/Dockerfile -t parsaghadermazi/myrnaseqpipeline:0.2.0 --push .
```

The default `docker` buildx driver **cannot** produce multi-platform manifests —
it fails with `Multi-platform build is not supported for the docker driver`.

Cell Ranger is deliberately absent — 10x's licence forbids redistributing it.
Supply it with `--cellranger_path`.

## Documentation

| | |
|---|---|
| [usage.md](docs/usage.md) | running each stage, parameters, output contract |
| [reference.md](docs/reference.md) | references and index building |
| [COVERAGE.md](docs/COVERAGE.md) | every tutorial tool vs what is implemented |
| [testing.md](docs/testing.md) | synthetic fixtures and what has been verified |
| [python-r-bridge.md](docs/python-r-bridge.md) | why no R module writes `.h5ad` |
| [design.md](design.md) | design principles and why |

### Lecture notes

The concepts behind each stage, written up as they were worked through:

- [01 — FASTQ to count matrix](docs/lectures/01-fastq-to-counts.md)
- [02 — QC and cell filtering](docs/lectures/02-qc.md)
- [03 — Integration](docs/lectures/03-integration.md)

## Testing

```bash
tests/run_tests.sh        # stage counts — 11/11 passing on v0.2.0
tests/run_guard_tests.sh  # refusals     — 2/2 passing
tests/run_qc_tests.sh     # stage qc     — 11/11 passing on v0.2.0
```

`run_guard_tests.sh` covers the failures that matter most: a design where batch
is perfectly confounded with condition, and a missing `.raw.h5ad`. Both must be
**refused** — a run that silently removes the treatment effect and emits a
clean-looking object is far more dangerous than a crash.

Both build synthetic fixtures with known ground truth, so results are checked by
equality rather than by eye.
