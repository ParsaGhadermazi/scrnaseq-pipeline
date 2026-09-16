#!/usr/bin/env python3
"""Normalisation and highly-variable gene selection on the merged object.

Library-size normalisation assumes a cell's total UMI count is technical. It is
partly biological -- plasma cells and proliferating cells genuinely hold more
RNA -- so this step converts expression into PROPORTIONS and absolute
transcriptional output is lost. That is a known, accepted cost of 'lognorm';
'pearson_residuals' avoids the naive division.

HVG selection is batch-aware. Selecting on pooled data picks genes that vary
BETWEEN samples -- donor HLA, sex-linked genes, dissociation-stress genes like
FOS/JUN/HSPA1A -- which then define PCA, the kNN graph, and therefore clusters.
Passing batch_key requires a gene to be variable WITHIN samples.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import scanpy as sc

__version__ = "0.1.0"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--method", default="lognorm", choices=["lognorm", "pearson_residuals"])
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--batch-key", default="sample")
    ap.add_argument("--n-pcs", type=int, default=50)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args()

    adata = sc.read_h5ad(a.input)
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()      # DE must always reach raw counts

    bk = a.batch_key if a.batch_key in adata.obs.columns else None
    if bk is None:
        print(f"WARNING: batch_key '{a.batch_key}' not in obs; HVG selection will NOT be batch-aware")

    if a.method == "pearson_residuals":
        sc.experimental.pp.highly_variable_genes(adata, flavor="pearson_residuals",
                                                 n_top_genes=a.n_hvg, batch_key=bk)
        sc.experimental.pp.normalize_pearson_residuals(adata)
    else:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=min(a.n_hvg, adata.n_vars - 1), batch_key=bk)

    adata.raw = adata
    sc.pp.scale(adata, max_value=10, zero_center=True)
    n_comps = int(min(a.n_pcs, adata.n_vars - 1, adata.n_obs - 1))
    sc.tl.pca(adata, n_comps=n_comps, mask_var="highly_variable", random_state=a.seed)

    adata.uns["normalisation"] = json.dumps({
        "method": a.method, "n_hvg": int(adata.var.highly_variable.sum()),
        "batch_key": bk, "n_pcs": n_comps})
    adata.write_h5ad(a.out)
    print(f"normalize_hvg: {a.method} | {adata.var.highly_variable.sum()} HVGs "
          f"(batch_key={bk}) | {n_comps} PCs")


if __name__ == "__main__":
    main()
