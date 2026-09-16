#!/usr/bin/env python3
"""Apply emptyDrops calls to the raw matrix, producing the cell-containing matrix.

Exists because R cannot write .h5ad (docs/python-r-bridge.md): emptydrops.R emits
per-barcode calls as CSV and this applies them.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import anndata as ad, pandas as pd

__version__ = "0.1.0"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw", required=True, type=Path)
    ap.add_argument("--calls", required=True, type=Path, help="ed_calls.csv from emptydrops.R")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--report", type=Path)
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args()

    adata = ad.read_h5ad(a.raw)
    calls = pd.read_csv(a.calls)
    keep = set(calls.loc[calls["is_cell"].astype(bool), "barcode"].astype(str))

    mask = adata.obs_names.astype(str).isin(keep)
    if mask.sum() == 0:
        raise SystemExit(f"ERROR: none of the {len(keep)} called cells matched barcodes in {a.raw}")

    out = adata[mask].copy()
    rep = {"n_barcodes_in": int(adata.n_obs), "n_cells_out": int(out.n_obs)}
    out.uns["emptydrops"] = json.dumps(rep)
    out.write_h5ad(a.out)
    print(f"subset_cells: {adata.n_obs} barcodes -> {out.n_obs} cells")

    if a.report:
        a.report.write_text(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()
