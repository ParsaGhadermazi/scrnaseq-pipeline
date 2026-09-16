#!/usr/bin/env python3
"""Preliminary, THROWAWAY clustering, used only to give SoupX its cluster labels.

This is the loop in Part 2's DAG: SoupX needs clusters, clustering needs
normalization, but the real normalization is supposed to happen AFTER ambient
correction. Resolution: cluster once on uncorrected data, use it only to label
cells for SoupX, then discard. Nothing downstream ever sees these labels.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import scanpy as sc, pandas as pd

__version__ = "0.1.0"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path, help="CSV: barcode,cluster")
    ap.add_argument("--resolution", type=float, default=0.8)
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--n-hvg", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args()

    adata = sc.read_h5ad(a.input)
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    n_hvg = min(a.n_hvg, adata.n_vars - 1)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg)
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=min(a.n_pcs, adata.n_vars - 1, adata.n_obs - 1), random_state=a.seed)
    sc.pp.neighbors(adata, n_pcs=min(a.n_pcs, adata.obsm["X_pca"].shape[1]), random_state=a.seed)
    sc.tl.leiden(adata, resolution=a.resolution, random_state=a.seed, flavor="igraph", n_iterations=2)

    pd.DataFrame({"barcode": adata.obs_names.astype(str),
                  "cluster": adata.obs["leiden"].astype(str)}).to_csv(a.out, index=False)
    print(f"prelim_cluster: {adata.n_obs} cells -> {adata.obs['leiden'].nunique()} clusters "
          f"(resolution {a.resolution}, discarded after SoupX)")


if __name__ == "__main__":
    main()
