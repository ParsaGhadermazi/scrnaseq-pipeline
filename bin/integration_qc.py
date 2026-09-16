#!/usr/bin/env python3
"""Score integration on two axes that pull against each other.

Each metric alone has a degenerate optimum: maximise mixing by collapsing every
cell into one blob, or maximise conservation by doing nothing at all. Only
reporting both makes over- and under-correction distinguishable.

  batch silhouette      LOW  is good (samples no longer separable)
  label silhouette      HIGH is good (cell populations still distinct)
  kNN batch mixing      HIGH is good (neighbourhoods contain mixed samples)

CAVEAT: true cell-type labels do not exist yet -- annotation is Part 4. Clusters
computed on the UNINTEGRATED embedding stand in for them. Those clusters are
themselves partly batch-driven, so the label silhouette catches gross
over-correction but is not ground truth. Stated in the output.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np, scanpy as sc
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

__version__ = "0.1.0"


def knn_mixing(emb: np.ndarray, labels: np.ndarray, k: int = 50, seed: int = 42) -> float:
    """Mean normalised entropy of batch labels among each cell's k neighbours.
    1.0 = neighbourhoods match global batch proportions; 0.0 = pure per-sample."""
    k = int(min(k, len(labels) - 1))
    nn = NearestNeighbors(n_neighbors=k + 1).fit(emb)
    idx = nn.kneighbors(emb, return_distance=False)[:, 1:]
    # np.unique(..., return_inverse=True) returns (values, inverse). Reversing
    # that tuple bound `uniq` to the VALUES (length = n_batches) and then indexed
    # it with neighbour indices up to n_obs:
    #   IndexError: index 8 is out of bounds for axis 0 with size 4
    values, codes = np.unique(labels, return_inverse=True)
    nb = codes[idx]                      # batch code of each neighbour
    counts = np.stack([(nb == i).sum(axis=1) for i in range(len(values))], axis=1)
    p = counts / counts.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -np.nansum(np.where(p > 0, p * np.log(p), 0.0), axis=1)
    return float(np.mean(ent) / np.log(len(values))) if len(values) > 1 else 1.0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--batch-key", default="sample")
    ap.add_argument("--n-neighbors", type=int, default=50)
    ap.add_argument("--subsample", type=int, default=5000, help="silhouette is O(n^2)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--version", action="version", version=__version__)
    a = ap.parse_args()

    adata = sc.read_h5ad(a.input)
    rng = np.random.default_rng(a.seed)

    # stand-in cell labels from the UNINTEGRATED embedding
    tmp = adata.copy()
    sc.pp.neighbors(tmp, use_rep="X_pca", random_state=a.seed)
    sc.tl.leiden(tmp, resolution=0.8, key_added="proxy", random_state=a.seed,
                 flavor="igraph", n_iterations=2)
    adata.obs["proxy_label"] = tmp.obs["proxy"].values

    res: dict = {"batch_key": a.batch_key,
                 "label_source": "leiden on UNINTEGRATED X_pca (proxy, not ground truth)",
                 "n_cells": int(adata.n_obs), "metrics": {}}

    idx = (rng.choice(adata.n_obs, a.subsample, replace=False)
           if adata.n_obs > a.subsample else np.arange(adata.n_obs))
    batch = adata.obs[a.batch_key].astype(str).to_numpy()
    label = adata.obs["proxy_label"].astype(str).to_numpy()

    for name, key in (("uncorrected", "X_pca"), ("integrated", "X_integrated")):
        if key not in adata.obsm:
            continue
        emb = np.asarray(adata.obsm[key])
        m = {}
        if len(np.unique(batch[idx])) > 1:
            m["batch_silhouette_lower_better"] = round(float(silhouette_score(emb[idx], batch[idx])), 4)
        if len(np.unique(label[idx])) > 1:
            m["label_silhouette_higher_better"] = round(float(silhouette_score(emb[idx], label[idx])), 4)
        m["knn_batch_mixing_higher_better"] = round(knn_mixing(emb, batch, a.n_neighbors, a.seed), 4)
        res["metrics"][name] = m

    u, i = res["metrics"].get("uncorrected", {}), res["metrics"].get("integrated", {})
    if u and i:
        res["delta"] = {
            "batch_silhouette": round(i.get("batch_silhouette_lower_better", 0)
                                      - u.get("batch_silhouette_lower_better", 0), 4),
            "label_silhouette": round(i.get("label_silhouette_higher_better", 0)
                                      - u.get("label_silhouette_higher_better", 0), 4),
            "knn_mixing": round(i.get("knn_batch_mixing_higher_better", 0)
                                - u.get("knn_batch_mixing_higher_better", 0), 4)}
        res["interpretation"] = (
            "Good integration: batch_silhouette DOWN, knn_mixing UP, label_silhouette roughly held. "
            "A large drop in label_silhouette alongside high mixing indicates OVER-correction -- "
            "populations merged, not aligned.")

    a.out.write_text(json.dumps(res, indent=2))
    print(json.dumps(res.get("delta", res["metrics"]), indent=2))


if __name__ == "__main__":
    main()
