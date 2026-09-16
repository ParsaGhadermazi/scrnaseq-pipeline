#!/usr/bin/env Rscript
# Ambient RNA ("soup") correction with SoupX.
#
# Needs BOTH matrices: the empty droplets ARE the soup measurement, so the raw
# matrix is not archival -- this step cannot run without it.
#
# Writes MTX + TSVs, never .h5ad (see docs/python-r-bridge.md).
suppressPackageStartupMessages({
  library(zellkonverter); library(SingleCellExperiment); library(SoupX); library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
get <- function(k, d = NULL) { i <- match(k, args); if (is.na(i)) d else args[i + 1] }

raw_h5ad  <- get("--raw")
cells_h5ad<- get("--cells")
clusters  <- get("--clusters")          # CSV: barcode,cluster  (preliminary, throwaway)
outdir    <- get("--outdir", "soupx_out")
rho_out   <- get("--rho", "rho.json")
force_rho <- get("--rho-force", "")      # override autoEstCont
min_rho   <- as.numeric(get("--min-rho", "0.05"))  # tutorial: correct only if >= 5%
tfidf_min <- as.numeric(get("--tfidf-min", "1.0"))
soup_q    <- as.numeric(get("--soup-quantile", "0.90"))
set.seed(as.integer(get("--seed", "42")))

grab <- function(f) {
  s <- readH5AD(f, reader = "R")
  assay(s, if ("counts" %in% assayNames(s)) "counts" else assayNames(s)[1])
}
tod <- grab(raw_h5ad)     # table of droplets -- ALL barcodes, incl. empties
toc <- grab(cells_h5ad)   # table of counts   -- real cells only

sc <- SoupChannel(tod, toc)

cl <- read.csv(clusters, stringsAsFactors = FALSE)
cl <- cl[match(colnames(toc), cl$barcode), ]
if (any(is.na(cl$cluster))) { cat("ERROR: clusters missing for some cells\n", file = stderr()); quit(status = 1) }
# Clusters are what make contamination IDENTIFIABLE: they let SoupX find genes a
# group of cells should NOT express, so any counts there are definitionally soup.
sc <- setClusters(sc, setNames(as.character(cl$cluster), cl$barcode))

if (nzchar(force_rho)) {
  rho <- as.numeric(force_rho); sc <- setContaminationFraction(sc, rho); est <- "manual"
} else {
  # autoEstCont needs genes a cluster should NOT express -- if every gene is
  # expressed everywhere ("low complexity") it finds no markers and stops.
  ok <- tryCatch({
    sc <- autoEstCont(sc, doPlot = FALSE, tfidfMin = tfidf_min, soupQuantile = soup_q); TRUE
  }, error = function(e) { cat("autoEstCont failed:", conditionMessage(e), "\n", file = stderr()); FALSE })
  if (!ok) {
    cat("\nSoupX could not estimate contamination automatically.\n",
        "  It needs marker genes: strongly expressed in some clusters, ABSENT from others.\n",
        "  Options: lower --tfidf-min (now ", tfidf_min, "), lower --soup-quantile (now ", soup_q, "),\n",
        "  or set --rho-force / params.soupx_rho from a known contamination estimate.\n",
        "  Refusing to guess a contamination fraction.\n", sep = "", file = stderr())
    quit(status = 1)
  }
  rho <- mean(sc$metaData$rho); est <- "autoEstCont"
}
cat("contamination rho = ", round(rho, 4), " (", est, ")\n", sep = "")

applied <- rho >= min_rho
out <- if (applied) adjustCounts(sc, roundToInt = TRUE) else toc
if (!applied) cat("rho < ", min_rho, " -- passing counts through uncorrected\n", sep = "")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
writeMM(out, file.path(outdir, "matrix.mtx"))
write.table(colnames(out), file.path(outdir, "barcodes.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rownames(out), file.path(outdir, "features.tsv"), quote = FALSE, row.names = FALSE, col.names = FALSE)

writeLines(sprintf('{"rho": %f, "estimator": "%s", "min_rho": %f, "applied": %s, "method": "soupx"}',
                   rho, est, min_rho, tolower(as.character(applied))), rho_out)
