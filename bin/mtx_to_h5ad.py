#!/usr/bin/env python3
"""Convert SoupX's MTX output back to .h5ad, restoring obs/var from the input.

The R side writes Matrix Market because no native R .h5ad writer exists
(docs/python-r-bridge.md). SoupX returns only the corrected counts, so cell and
gene metadata are carried over from the pre-correction object rather than lost.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp

__version__ = "0.1.0"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mtx-dir", required=True, type=Path)
    ap.add_argument("--source", required=True, type=Path, help="pre-correction .h5ad, for obs/var")
    ap.add_argument("--rho", type=Path, help="rho.json from SoupX")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args()

    m = sio.mmread(a.mtx_dir / "matrix.mtx").tocsr()            # genes x cells
    barcodes = [l.strip() for l in open(a.mtx_dir / "barcodes.tsv") if l.strip()]
    features = [l.strip() for l in open(a.mtx_dir / "features.tsv") if l.strip()]

    src = ad.read_h5ad(a.source)
    out = ad.AnnData(m.T.tocsr(),                                # -> cells x genes
                     obs=src.obs.reindex(barcodes).copy(),
                     var=src.var.reindex(features).copy())
    out.obs_names, out.var_names = barcodes, features

    # keep the uncorrected counts: DE (Part 5) must never run on corrected values
    pre = src[barcodes, features]
    out.layers["counts_raw"] = pre.X.copy()
    out.X = out.X.astype(np.int32)

    # Copy only values h5py can actually write. Upstream uns entries are JSON
    # strings by convention (see cell_qc.py); anything else is stringified rather
    # than risking a write failure three processes later.
    for k, v in src.uns.items():
        out.uns[k] = v if isinstance(v, (str, int, float, bool)) else json.dumps(v, default=str)
    rho_info = {}
    if a.rho and a.rho.exists():
        raw_rho = a.rho.read_text().strip()
        out.uns["ambient"] = raw_rho          # string: h5py cannot write a dict
        try:
            rho_info = json.loads(raw_rho)    # dict: for the log line only
        except json.JSONDecodeError:
            rho_info = {}

    out.write_h5ad(a.out)
    print(f"mtx_to_h5ad: {out.n_obs} cells x {out.n_vars} genes | "
          f"rho={rho_info.get('rho', 'n/a')} applied={rho_info.get('applied', 'n/a')}")


if __name__ == "__main__":
    main()
