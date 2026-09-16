#!/usr/bin/env Rscript
# Distinguish real cells from empty droplets on the RAW matrix.
#
# Emits a CSV of per-barcode calls, NOT a matrix: R cannot write .h5ad
# (see docs/python-r-bridge.md). subset_cells.py applies these calls.
suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(DropletUtils)
})

args <- commandArgs(trailingOnly = TRUE)
get <- function(k, d = NULL) { i <- match(k, args); if (is.na(i)) d else args[i + 1] }

raw_h5ad <- get("--raw")
out_csv  <- get("--out")
lower    <- as.integer(get("--lower", "100"))
niters   <- as.integer(get("--niters", "10000"))
fdr_max  <- as.numeric(get("--fdr", "0.01"))
seed     <- as.integer(get("--seed", "42"))
set.seed(seed)                       # emptyDrops is Monte Carlo: seed it or results move

# reader="R" is essential -- the default goes through basilisk and compiles
# Python at runtime. See docs/python-r-bridge.md.
sce <- readH5AD(raw_h5ad, reader = "R")
assay_name <- if ("counts" %in% assayNames(sce)) "counts" else assayNames(sce)[1]
m <- assay(sce, assay_name)

cat("emptyDrops: ", nrow(m), " genes x ", ncol(m), " barcodes (lower=", lower,
    ", niters=", niters, ")\n", sep = "")

# Barcodes below `lower` are ASSUMED empty and become the ambient profile.
# Everything above is tested against it, so a low-count barcode survives if its
# profile is DISTINCTIVE -- which is what saves small cells (platelets,
# neutrophils) that a plain UMI-count knee would discard.
e <- emptyDrops(m, lower = lower, niters = niters)

res <- data.frame(
  barcode = colnames(m),
  total   = e$Total,
  pvalue  = e$PValue,
  fdr     = e$FDR,
  is_cell = !is.na(e$FDR) & e$FDR <= fdr_max
)
write.csv(res, out_csv, row.names = FALSE, quote = FALSE)

n <- sum(res$is_cell)
cat("cells called: ", n, " / ", nrow(res), " barcodes (FDR <= ", fdr_max, ")\n", sep = "")
if (n == 0) { cat("ERROR: emptyDrops called zero cells\n", file = stderr()); quit(status = 1) }
