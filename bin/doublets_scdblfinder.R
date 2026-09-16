#!/usr/bin/env Rscript
# Doublet detection with scDblFinder.
#
# Emits per-cell calls as CSV -- no matrix crosses the R/Python boundary.
# scDblFinder self-estimates the doublet rate, so unlike DoubletFinder it needs
# neither nExp nor a pK sweep.
suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(scDblFinder)
})

args <- commandArgs(trailingOnly = TRUE)
get <- function(k, d = NULL) { i <- match(k, args); if (is.na(i)) d else args[i + 1] }

in_h5ad <- get("--input")
out_csv <- get("--out")
clusters<- get("--clusters", "")
dbr     <- get("--dbr", "")        # blank => estimate from cell count (~0.8%/1000)
set.seed(as.integer(get("--seed", "42")))

sce <- readH5AD(in_h5ad, reader = "R")
if (!"counts" %in% assayNames(sce)) assayNames(sce)[1] <- "counts"

a <- list(sce = sce)
if (nzchar(clusters)) {
  cl <- read.csv(clusters, stringsAsFactors = FALSE)
  cl <- cl[match(colnames(sce), cl$barcode), ]
  a$clusters <- as.character(cl$cluster)
}
if (nzchar(dbr)) a$dbr <- as.numeric(dbr)

sce <- do.call(scDblFinder, a)

res <- data.frame(
  barcode        = colnames(sce),
  doublet_score  = sce$scDblFinder.score,
  doublet_class  = as.character(sce$scDblFinder.class),
  method         = "scDblFinder"
)
write.csv(res, out_csv, row.names = FALSE, quote = FALSE)

# NOTE: this is mostly blind to HOMOTYPIC doublets (two cells of the same type),
# which look like one larger normal cell. That residual is carried downstream.
cat("doublets: ", sum(res$doublet_class == "doublet"), " / ", nrow(res),
    " (", round(100 * mean(res$doublet_class == "doublet"), 1), "%)\n", sep = "")
