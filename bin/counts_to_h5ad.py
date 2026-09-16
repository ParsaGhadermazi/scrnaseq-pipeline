#!/usr/bin/env python3
"""Normalise any quantifier's native output into one .h5ad schema.

This is the Part 1 / Part 2 boundary. Whichever quantifier ran, everything
downstream sees exactly this:

    X                     csr, cells x genes, RAW counts, never normalised
    obs.index             barcode, suffix stripped
    obs[...]              samplesheet passthrough (sample, condition, batch, ...)
    var.index             gene_id, Ensembl, version suffix stripped
    var['gene_symbol']    NOT unique -- never used as a key
    var['feature_type']   'Gene Expression' | 'Antibody Capture' | ...
    uns['pipeline']       quantifier, chemistry, rounding policy, tool versions

Six places quantifiers disagree, all absorbed here (see lecture 1, section 10):
  1. orientation  -- Cell Ranger MTX is genes x cells; AnnData is cells x genes
  2. barcode suffix -- Cell Ranger appends '-1'; alevin-fry and kb do not
  3. gene id version -- ENSG00000141510.17 vs ENSG00000141510
  4. gene universe  -- differs by reference; recorded in uns, never silently intersected
  5. cell calling   -- each tool's "filtered" means something slightly different
  6. FRACTIONAL counts -- alevin-fry's EM UMI resolution splits ambiguous UMIs
     across genes. DESeq2 and scDblFinder assume integers. The rounding policy
     is explicit and recorded, never a silent astype(int).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

__version__ = "0.1.0"


# --------------------------------------------------------------- readers
def read_cellranger(quant: Path):
    """quant = the `outs/` directory."""
    import scanpy as sc

    filt = quant / "filtered_feature_bc_matrix.h5"
    raw = quant / "raw_feature_bc_matrix.h5"
    if not filt.exists():
        sys.exit(f"ERROR: {filt} not found -- is --quant pointing at cellranger's outs/ ?")
    a = sc.read_10x_h5(filt, gex_only=False)
    b = sc.read_10x_h5(raw, gex_only=False) if raw.exists() else None
    return a, b


def read_alevinfry(quant: Path):
    """quant = simpleaf's output dir, containing af_quant/."""
    from pyroe import load_fry

    inner = quant / "af_quant" if (quant / "af_quant").exists() else quant
    # 'raw' keeps S/U/A layers separate so Part 13 (scVelo) stays possible;
    # 'snRNA' collapses them. We keep both: X = collapsed, layers = components.
    a = load_fry(str(inner), output_format="raw")
    if {"spliced", "unspliced", "ambiguous"} <= set(a.layers):
        a.X = a.layers["spliced"] + a.layers["ambiguous"]
    # --unfiltered-pl gives one permit-list matrix: it IS the raw matrix.
    return a, a.copy()


def read_kallisto(quant: Path):
    """quant = kb's output dir, with counts_filtered/ and counts_unfiltered/."""
    import scanpy as sc

    def _one(d: Path):
        if not d.exists():
            return None
        a = sc.read_mtx(d / "cells_x_genes.mtx")          # already cells x genes
        a.obs_names = (d / "cells_x_genes.barcodes.txt").read_text().split()
        a.var_names = (d / "cells_x_genes.genes.txt").read_text().split()
        # kb writes symbols to a separate file; without this every gene_symbol
        # would just be a copy of the gene_id.
        names = d / "cells_x_genes.genes.names.txt"
        if names.exists():
            sym = names.read_text().split()
            if len(sym) == a.n_vars:
                a.var["gene_symbol"] = sym
        return a

    return _one(quant / "counts_filtered"), _one(quant / "counts_unfiltered")


READERS = {"cellranger": read_cellranger, "alevinfry": read_alevinfry, "kallisto": read_kallisto}


# --------------------------------------------------------------- normalisation
def normalise(adata, obs_extra: dict, quantifier: str, rounding: str):
    import anndata  # noqa: F401

    # 2. barcode suffix
    adata.obs_names = [b.split("-")[0] for b in adata.obs_names]

    # 3. gene id version  +  keep symbols in a column, never as the key
    if "gene_ids" in adata.var.columns:                    # cellranger: var_names are symbols
        adata.var["gene_symbol"] = adata.var_names.astype(str)
        adata.var_names = adata.var["gene_ids"].astype(str)
        adata.var.drop(columns=["gene_ids"], inplace=True)
    elif "gene_symbol" not in adata.var.columns:
        adata.var["gene_symbol"] = adata.var_names.astype(str)
    adata.var_names = [g.split(".")[0] for g in adata.var_names]
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    if "feature_types" in adata.var.columns:
        adata.var["feature_type"] = adata.var.pop("feature_types").astype(str)
    if "feature_type" not in adata.var.columns:
        adata.var["feature_type"] = "Gene Expression"

    # 6. fractional counts
    X = adata.X
    was_fractional = bool(np.issubdtype(X.dtype, np.floating)
                          and not np.allclose(X.data if sp.issparse(X) else X,
                                              np.rint(X.data if sp.issparse(X) else X)))
    if was_fractional and rounding != "keep":
        fn = np.rint if rounding == "round" else np.floor
        if sp.issparse(X):
            X.data = fn(X.data)
            X.eliminate_zeros()
        else:
            X = fn(X)
    if rounding != "keep":
        X = X.astype(np.int32)
    adata.X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)

    # Readers add their own bookkeeping columns (pyroe.load_fry contributes
    # obs['barcodes'], a duplicate of the index). Drop anything that merely
    # restates obs_names, or the "identical schema across quantifiers"
    # guarantee quietly stops holding.
    for col in list(adata.obs.columns):
        if adata.obs[col].astype(str).equals(adata.obs_names.to_series().astype(str)):
            adata.obs.drop(columns=[col], inplace=True)

    for k, v in obs_extra.items():
        adata.obs[k] = v

    adata.uns["pipeline"] = {
        "quantifier": quantifier,
        "counts_to_h5ad_version": __version__,
        "rounding_policy": rounding,
        "counts_were_fractional": was_fractional,
        "gene_key": "ensembl_gene_id (version stripped)",
        **{k: v for k, v in obs_extra.items() if k in ("chemistry", "reference")},
    }
    return adata


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quant", required=True, type=Path)
    ap.add_argument("--quantifier", required=True, choices=sorted(READERS))
    ap.add_argument("--sample", required=True)
    ap.add_argument("--obs", default="{}", help="JSON of samplesheet metadata -> obs columns")
    ap.add_argument("--rounding", default="round", choices=["round", "floor", "keep"])
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--out-raw", type=Path)
    ap.add_argument("--version", action="version", version=__version__)
    args = ap.parse_args()

    obs_extra = json.loads(args.obs)
    obs_extra.setdefault("sample", args.sample)

    filtered, raw = READERS[args.quantifier](args.quant)
    if filtered is None:
        sys.exit(f"ERROR: no filtered matrix found under {args.quant}")

    normalise(filtered, obs_extra, args.quantifier, args.rounding).write_h5ad(args.out)
    print(f"wrote {args.out}  ({filtered.n_obs} cells x {filtered.n_vars} genes)")

    # The raw matrix is not optional bookkeeping: SoupX and emptyDrops in Part 2
    # profile the ambient soup from the EMPTY droplets, which exist only here.
    if args.out_raw and raw is not None:
        normalise(raw, obs_extra, args.quantifier, args.rounding).write_h5ad(args.out_raw)
        print(f"wrote {args.out_raw}  ({raw.n_obs} barcodes x {raw.n_vars} genes)")
    elif args.out_raw:
        print("WARNING: no raw matrix available -- SoupX/emptyDrops in Part 2 will be limited",
              file=sys.stderr)


if __name__ == "__main__":
    main()
