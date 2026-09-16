#!/usr/bin/env python3
"""Cell- and gene-level QC filtering.

Thresholds are either FIXED (reproduces the tutorial exactly) or MAD-based
(default; generalises across tissues and sequencing depths).

Why MAD and not mean/SD: standard deviation squares deviations, so the outlier
cells you are hunting inflate the very statistic used to find them, widening the
gate until they fall inside it. Median and MAD have a 50% breakdown point.

Gene symbols are matched against var['gene_symbol'] because var_names are
Ensembl IDs -- symbols are NOT unique and are never used as a key.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scanpy as sc

__version__ = "0.1.0"


def symbols(adata) -> pd.Series:
    col = "gene_symbol" if "gene_symbol" in adata.var.columns else None
    return (adata.var[col] if col else pd.Series(adata.var_names, index=adata.var_names)).astype(str)


def mad_bounds(x: np.ndarray, nmads: float, log: bool) -> tuple[float, float, bool]:
    """Return (lo, hi, degenerate). MAD==0 when >50% of values are identical --
    real in sparse data (percent_mt is often exactly 0 in nuclei). Dividing by it
    would flag every non-identical cell, so we report it instead of filtering."""
    v = np.log1p(x) if log else x
    med = float(np.median(v))
    mad = float(np.median(np.abs(v - med)))
    if mad == 0:
        return (-np.inf, np.inf, True)
    mad *= 1.4826                                   # scale to be comparable to SD
    lo, hi = med - nmads * mad, med + nmads * mad
    return (float(np.expm1(lo)) if log else lo, float(np.expm1(hi)) if log else hi, False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--doublets", type=Path, nargs="*", default=[], help="doublet call CSVs")
    ap.add_argument("--doublet-removal", default="primary", choices=["primary", "consensus", "none"])
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--report", required=True, type=Path)
    ap.add_argument("--method", default="mad", choices=["mad", "fixed"])
    ap.add_argument("--nmads", type=float, default=3.0)
    ap.add_argument("--min-genes", type=int, default=500)
    ap.add_argument("--max-genes", type=int, default=5000)
    ap.add_argument("--min-counts", type=int, default=800)
    ap.add_argument("--max-counts", type=int, default=20000)
    ap.add_argument("--max-mito-pct", type=float, default=20.0, help="absolute ceiling, always applied")
    ap.add_argument("--min-cells-per-gene", type=int, default=3)
    ap.add_argument("--remove-hb-genes", default="true")
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args()

    adata = ad.read_h5ad(a.input)
    n0 = adata.n_obs
    rep: dict = {"method": a.method, "n_cells_in": int(n0), "n_genes_in": int(adata.n_vars), "steps": []}

    sym = symbols(adata)
    adata.var["mt"]   = sym.str.upper().str.startswith("MT-").values
    adata.var["ribo"] = sym.str.upper().str.match(r"^RP[SL]").fillna(False).values
    adata.var["hb"]   = sym.str.upper().str.match(r"^HB[ABDEGMQZ]\d*$").fillna(False).values
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=False, percent_top=None)

    # ---- doublets -----------------------------------------------------------
    if a.doublets and a.doublet_removal != "none":
        calls = [pd.read_csv(f) for f in a.doublets]
        methods = [c["method"].iloc[0] if "method" in c else f"m{i}" for i, c in enumerate(calls)]
        wide = pd.DataFrame(index=adata.obs_names.astype(str))
        for m, c in zip(methods, calls):
            wide[m] = c.set_index(c["barcode"].astype(str))["doublet_class"].reindex(wide.index).eq("doublet")
        adata.obs["doublet_n_methods"] = wide.sum(axis=1).values
        is_dbl = wide.all(axis=1) if a.doublet_removal == "consensus" else wide[methods[0]]
        adata.obs["is_doublet"] = is_dbl.fillna(False).values
        before = adata.n_obs
        adata = adata[~adata.obs["is_doublet"]].copy()
        rep["steps"].append({"step": "doublets", "rule": a.doublet_removal, "methods": methods,
                             "removed": int(before - adata.n_obs), "remaining": int(adata.n_obs)})

    # ---- cell QC ------------------------------------------------------------
    degenerate = []
    if a.method == "mad":
        lo_c, hi_c, d1 = mad_bounds(adata.obs["total_counts"].to_numpy(), a.nmads, log=True)
        lo_g, hi_g, d2 = mad_bounds(adata.obs["n_genes_by_counts"].to_numpy(), a.nmads, log=True)
        _,   hi_m, d3 = mad_bounds(adata.obs["pct_counts_mt"].to_numpy(), a.nmads, log=False)
        for name, d in (("total_counts", d1), ("n_genes", d2), ("pct_mt", d3)):
            if d:
                degenerate.append(name)
        # the absolute mito ceiling always applies: in a very clean sample MAD is
        # tiny and would cut healthy cells; in a dying one the median itself is high
        hi_m = min(hi_m, a.max_mito_pct) if np.isfinite(hi_m) else a.max_mito_pct
    else:
        lo_c, hi_c, lo_g, hi_g, hi_m = a.min_counts, a.max_counts, a.min_genes, a.max_genes, a.max_mito_pct

    keep = ((adata.obs["total_counts"] >= lo_c) & (adata.obs["total_counts"] <= hi_c) &
            (adata.obs["n_genes_by_counts"] >= lo_g) & (adata.obs["n_genes_by_counts"] <= hi_g) &
            (adata.obs["pct_counts_mt"] <= hi_m))
    before = adata.n_obs
    adata = adata[keep].copy()
    rep["thresholds"] = {"total_counts": [lo_c, hi_c], "n_genes": [lo_g, hi_g], "pct_mt_max": hi_m,
                         "degenerate_mad": degenerate}
    rep["steps"].append({"step": "cell_qc", "removed": int(before - adata.n_obs), "remaining": int(adata.n_obs)})
    if degenerate:
        print(f"WARNING: MAD was 0 for {degenerate} (>50% identical values); "
              f"those metrics were NOT filtered on", flush=True)

    # ---- gene QC ------------------------------------------------------------
    g0 = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=a.min_cells_per_gene)
    if str(a.remove_hb_genes).lower() == "true":
        # RBC contamination in PBMCs -- but in bone marrow / erythropoiesis these
        # genes ARE the biology, which is why this is a parameter.
        adata = adata[:, ~adata.var["hb"].astype(bool)].copy()
    rep["steps"].append({"step": "gene_qc", "removed": int(g0 - adata.n_vars), "remaining": int(adata.n_vars)})

    rep.update(n_cells_out=int(adata.n_obs), n_genes_out=int(adata.n_vars),
               pct_cells_kept=round(100 * adata.n_obs / max(n0, 1), 2))
    # uns must be h5ad-writable: rep contains a list of dicts, which h5py cannot
    # serialise ("Can't implicitly convert non-string objects to strings").
    # Keep the machine-readable form in the sidecar json; store a string here.
    adata.uns["cell_qc"] = json.dumps(rep)
    adata.write_h5ad(a.out)
    a.report.write_text(json.dumps(rep, indent=2))
    print(f"cell_qc: {n0} -> {adata.n_obs} cells ({rep['pct_cells_kept']}%), {g0} -> {adata.n_vars} genes")


if __name__ == "__main__":
    main()
