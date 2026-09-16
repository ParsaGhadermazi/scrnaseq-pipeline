#!/usr/bin/env python3
"""Batch integration.

Integration corrects an EMBEDDING, not expression. Counts and log-normalised
values are left untouched, so Part 5's differential expression reads
layers['counts'] and never sees corrected values -- feeding shifted continuous
values to a negative-binomial model breaks it, and integration deliberately
removes the between-sample variance DE needs to estimate uncertainty.

Two guards run before any correction:
  * one batch level        -> nothing to integrate; skip rather than fabricate
  * batch == condition     -> refuse; correcting would delete the biology
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np, scanpy as sc, pandas as pd

__version__ = "0.1.0"
R_ONLY = {"cca", "rpca", "fastmnn"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--method", default="harmony")
    ap.add_argument("--batch-key", default="sample")
    ap.add_argument("--condition-key", default="condition")
    ap.add_argument("--n-pcs", type=int, default=30)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args()

    adata = sc.read_h5ad(a.input)
    m = a.method.lower()
    if m in R_ONLY:
        sys.exit(f"ERROR: '{m}' is a Seurat method and runs in the R module, not here")

    if a.batch_key not in adata.obs.columns:
        sys.exit(f"ERROR: batch_key '{a.batch_key}' not found in obs")

    nb = adata.obs[a.batch_key].nunique()
    if m == "none" or nb < 2:
        adata.obsm["X_integrated"] = adata.obsm["X_pca"].copy()
        adata.uns["integration"] = json.dumps({
            "method": "none",
            "reason": f"only {nb} batch level(s)" if nb < 2 else "requested",
            "batch_key": a.batch_key, "embedding": "X_integrated"})
        adata.write_h5ad(a.out)
        print(f"integrate: skipped ({nb} batch level(s)); X_integrated = X_pca")
        return

    # confounding guard: if every batch maps to exactly one condition and vice
    # versa, the two variables are the same and correction removes the result
    if a.condition_key in adata.obs.columns:
        ct = pd.crosstab(adata.obs[a.batch_key], adata.obs[a.condition_key])
        if (ct > 0).sum(axis=1).max() == 1 and adata.obs[a.condition_key].nunique() > 1:
            per_cond = (ct > 0).sum(axis=0)
            if (per_cond == 1).all():
                sys.exit(
                    f"ERROR: '{a.batch_key}' is perfectly confounded with '{a.condition_key}'.\n"
                    "  Each batch contains exactly one condition and each condition one batch, so\n"
                    "  they are the same variable. Integrating would remove the condition effect.\n"
                    "  This is an experimental design problem; no setting fixes it.")

    if m == "harmony":
        # Call harmonypy directly rather than sc.external.pp.harmony_integrate.
        # That wrapper does `harmony_out.Z_corr.T`, which was right for
        # harmonypy 1.x (Z_corr was PCs x cells). harmonypy 2.0.0 returns
        # Z_corr as cells x PCs, so the wrapper transposes it to (n_pcs, n_obs)
        # and AnnData rejects it:
        #   "Value passed for key 'X_integrated' is of incorrect shape...
        #    Value had shape (50,) while it should have had (1930,)"
        # Orient on n_obs instead of trusting either convention.
        import harmonypy
        X = np.asarray(adata.obsm["X_pca"], dtype=np.float64)
        ho = harmonypy.run_harmony(X, adata.obs, a.batch_key, random_state=a.seed)
        Z = np.asarray(ho.Z_corr)
        if Z.shape[0] == adata.n_obs:
            emb = Z
        elif Z.shape[1] == adata.n_obs:
            emb = Z.T
        else:
            sys.exit(f"ERROR: harmony returned {Z.shape}; neither axis matches "
                     f"n_obs={adata.n_obs}")
        adata.obsm["X_integrated"] = emb
    elif m == "scanorama":
        sc.external.pp.scanorama_integrate(adata, key=a.batch_key, basis="X_pca",
                                           adjusted_basis="X_integrated")
    else:
        sys.exit(f"ERROR: unknown integration method '{m}'")

    adata.uns["integration"] = json.dumps({
        "method": m, "batch_key": a.batch_key, "n_batches": int(nb),
        "embedding": "X_integrated",
        "note": "embedding only; counts and X are unmodified"})
    adata.write_h5ad(a.out)
    print(f"integrate: {m} over {nb} '{a.batch_key}' levels -> obsm['X_integrated']")


if __name__ == "__main__":
    main()
