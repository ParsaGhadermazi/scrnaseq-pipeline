# Single image for the whole pipeline, built multi-arch:
#   docker buildx build --platform linux/amd64,linux/arm64 \
#     -f containers/Dockerfile -t parsaghadermazi/myrnaseqpipeline:0.2.0 --push .
#
# The tag MUST match params.container_image in nextflow.config. Nextflow pulls by
# that exact name, so a mismatch fails at container launch -- and `-preview`
# cannot catch it, because it never starts a container.
#
# Grown one layer per lecture:
#   v0.1.0  Part 1  -- fastq -> counts, anndata core
#   v0.2.0  Parts 2/3 -- + R stack (Seurat, SoupX, scDblFinder, harmony, ...)
#
# cellranger is NOT here: 10x's licence forbids redistribution. It is
# bind-mounted at runtime from params.cellranger_path.
FROM mambaorg/micromamba:1.5-noble

ARG TARGETARCH
LABEL org.opencontainers.image.title="myRNASeqPipeline" \
      org.opencontainers.image.version="0.2.0" \
      org.opencontainers.image.source="https://hub.docker.com/r/parsaghadermazi/myrnaseqpipeline" \
      org.opencontainers.image.description="scRNA-seq pipeline: fastq to counts, QC, integration"

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git curl ca-certificates pkg-config \
        zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
        libboost-all-dev libhdf5-dev autoconf flex bison \
    && rm -rf /var/lib/apt/lists/*

ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=/opt/conda/bin:$PATH

# 1. arch-safe core
COPY env/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -afy

# 2. arch-fragile quantifiers, conda-first with source fallback
COPY containers/install_quantifiers.sh /tmp/install_quantifiers.sh
COPY env/environment.quant.yml /tmp/environment.quant.yml
RUN bash /tmp/install_quantifiers.sh /tmp/environment.quant.yml

# 3. R stack for Part 2 (QC) and Part 3 (integration).
#    THE .h5ad BRIDGE IS ASYMMETRIC -- verified, not assumed:
#      READ  : readH5AD(file, reader = "R")  works offline via rhdf5.  OK
#      WRITE : there is NO native R writer. `writer` is not an argument of
#              writeH5AD at all (formals: sce, file, X_name, skip_assays,
#              compression, version, ...). Passing it is silently swallowed by
#              `...` and the PYTHON path runs, which bootstraps basilisk/pyenv
#              and compiles Python at RUNTIME. Non-reproducible, fails offline.
#
#    Therefore: R modules READ .h5ad directly, and any R step that must return a
#    MATRIX writes Matrix::writeMM (MTX) + TSVs, which a small Python step
#    converts to .h5ad. Only ambient correction rewrites a matrix; every other
#    R step emits per-cell CSVs, so this costs one conversion in the whole DAG.
COPY containers/install_r_packages.sh /tmp/install_r_packages.sh
COPY env/environment.r.yml /tmp/environment.r.yml
RUN bash /tmp/install_r_packages.sh /tmp/environment.r.yml

# 4. pipeline scripts are NOT baked in.
#    Nextflow stages $projectDir/bin onto every task's PATH automatically, so
#    baking them here would only create a stale-code trap: the image copy would
#    shadow the repo copy and silently run an old script. Script versions are
#    pinned by the pipeline's own git revision instead.
RUN mkdir -p /opt/pipeline

# 5. freeze what actually resolved -- this file, not environment.yml,
#    is the authoritative version record.
RUN micromamba env export -n base > /opt/pipeline/environment.lock.yml

WORKDIR /work
