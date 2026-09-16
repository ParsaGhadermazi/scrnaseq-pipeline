#!/usr/bin/env python3
"""Merge per-sample QC'd .h5ad files into one object.

Genes are intersected, not unioned: a union would fabricate structured zeros for
genes absent from a sample's reference, and those zeros look like biology to
HVG selection and PCA.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import anndata as ad

__version__ = "0.1.0"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", required=True, nargs="+", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--report", type=Path)
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args()

    adatas = {}
    for f in a.inputs:
        x = ad.read_h5ad(f)
        name = str(x.obs["sample"].iloc[0]) if "sample" in x.obs.columns else f.stem
        adatas[name] = x

    merged = ad.concat(adatas, join="inner", label="sample_key", index_unique="-", merge="same")
    genes = {k: v.n_vars for k, v in adatas.items()}
    rep = {"n_samples": len(adatas), "cells_per_sample": {k: int(v.n_obs) for k, v in adatas.items()},
           "genes_per_sample": {k: int(v) for k, v in genes.items()},
           "genes_after_intersection": int(merged.n_vars), "n_cells": int(merged.n_obs)}
    # same h5ad constraint as cell_qc.py: rep holds nested dicts
    merged.uns["merge"] = json.dumps(rep)
    merged.write_h5ad(a.out)
    if a.report:
        a.report.write_text(json.dumps(rep, indent=2))
    print(f"merge: {len(adatas)} samples -> {merged.n_obs} cells x {merged.n_vars} genes "
          f"(intersection; per-sample range {min(genes.values())}-{max(genes.values())})")


if __name__ == "__main__":
    main()
