#!/usr/bin/env bash
# Build a synthetic fixture for stage `qc`, with known ground truth.
#
# Stage qc needs more structure than stage counts: emptyDrops needs EMPTY
# droplets to profile the soup, SoupX needs contamination actually present,
# scDblFinder needs doublets to find, and integration needs a batch effect to
# remove. All of that is planted deliberately so results can be checked by
# equality rather than by eye.
#
#   usage: tests/make_qc_fixture.sh <outdir> [image]
set -euo pipefail
OUT="${1:?usage: make_qc_fixture.sh <outdir> [image]}"
IMG="${2:-parsaghadermazi/myrnaseqpipeline:0.2.0}"

if ! docker image inspect "$IMG" >/dev/null 2>&1; then
  echo "ERROR: image '$IMG' not found locally (built, never pushed)." >&2
  exit 1
fi
mkdir -p "$OUT/h5ad"; OUT="$(cd "$OUT" && pwd)"

docker run --rm -i -v "$OUT":/w -w /w "$IMG" python - <<'PY'
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp, pathlib, json
rng = np.random.default_rng(7)

# Gene count and ambient fraction are not cosmetic. SoupX's quickMarkers ranks
# by tf-idf, and idf = log(n_cells / n_cells_expressing_gene). With few genes,
# ambient deposits >1 count of EVERY gene in EVERY cell, every gene is detected
# everywhere, idf collapses to ~0, and NO gene passes tfidfMin -- however
# exclusive its biological expression is. A 400-gene fixture gave markers
# idf=0.083 vs housekeeping idf=0.083: zero discriminative power.
# Real data has ~30k genes, so per-gene ambient is <<1 count and off-target
# markers are genuinely zero in most cells. The fixture must match that SPARSITY,
# not merely the fold-change.
#
# THE BINDING CONSTRAINT IS N_TYPES, measured by sweep (not derived):
#
#     idf <= log(N_TYPES)
#
# A PERFECTLY exclusive marker is still expressed in 1/N_TYPES of all cells, so
# idf can never exceed log(N_TYPES) even with zero ambient. Verified: at ambient=0
# measured idf equals log(N_TYPES) to 3 decimals for N_TYPES in {3,6,8,12}.
#
# With 3 types the ceiling is log(3) = 1.099, and ANY ambient pushes idf below
# SoupX's tfidfMin of 1.0. No value of MARKER_LEVEL or AMBIENT_FRAC can rescue a
# 3-type fixture -- which is exactly why two rounds of tuning those moved idf
# only 0.57 -> 0.66.
#
# Measured idf for type-0 markers:
#     n_types   ambient=0.08   ambient=0.03
#         3         0.728          0.930
#         6         1.313          1.575
#         8         1.585          1.845
#        12         1.948          2.242
#
# 12 types / MARKER_LEVEL 30 gives idf ~1.38, clearing tfidfMin=1.0 with margin.
#
# The SECOND gate (soupQuantile) is NOT reachable from this fixture. A joint
# sweep over n_types {8,12,16} x marker_level {30,40} x hk_scale {0.25,0.10}
# left soup/q90 pinned at 0.98-0.99 in EVERY cell -- flat across all three
# levers, i.e. not tunable by them. With ~1300 genes a synthetic fixture cannot
# reproduce the soup-mass distribution of 30k-gene real data, where a marker of
# an abundant cell type genuinely is high in the soup.
# So conf/test.config lowers soupx_soup_quantile for the TEST PROFILE ONLY;
# production keeps SoupX's default 0.90.
N_TYPES, MARKERS_PER_TYPE, HOUSEKEEPING = 12, 40, 980
N_MARKERS = N_TYPES * MARKERS_PER_TYPE
N_GENES   = N_MARKERS + HOUSEKEEPING
CELLS_PER_TYPE, N_EMPTY, DOUBLET_N, AMBIENT_FRAC = 40, 1000, 25, 0.08
DEPTH = 1.0

# CRITICAL: marker genes must be EXCLUSIVE -- expressed in their own cell type
# and ~absent elsewhere. SoupX's autoEstCont estimates contamination from genes a
# cluster should NOT express; if every gene is nonzero in every type there are no
# "should be zero" groups and it fails with "No plausible marker genes found".
# Real data has this property (CD3D in T cells, not monocytes); a fixture must too.
blocks = {t: np.arange(t*MARKERS_PER_TYPE, (t+1)*MARKERS_PER_TYPE) for t in range(N_TYPES)}
hk_slice = slice(N_MARKERS, N_GENES)
base_hk  = rng.gamma(2.0, 0.25, HOUSEKEEPING)    # expressed in every cell type
MARKER_LEVEL = 30.0

def profile(t, scale):
    mu = np.zeros(N_GENES)
    mu[hk_slice] = base_hk                        # housekeeping everywhere
    mu[blocks[t]] = MARKER_LEVEL                  # this type's markers only
    return mu * scale

def make_sample(name, batch_scale, seed, cond):
    r = np.random.default_rng(seed)
    cells, labels = [], []
    for t in range(N_TYPES):
        mu = profile(t, batch_scale)
        for _ in range(CELLS_PER_TYPE):
            cells.append(r.poisson(mu * DEPTH)); labels.append(f"type{t}")
    for _ in range(DOUBLET_N):                    # heterotypic: two marker blocks lit
        a, b = r.choice(N_TYPES, 2, replace=False)
        cells.append(r.poisson((profile(a, batch_scale) + profile(b, batch_scale)) * DEPTH * 0.5))
        labels.append("doublet")

    X = np.vstack(cells).astype(np.int32)
    # the soup is the average of all cells, so it CONTAINS marker genes -- which
    # is what puts type-0 markers into type-1 cells for SoupX to measure
    soup = X.sum(0).astype(float); soup /= soup.sum()
    X = X + r.poisson(np.outer(X.sum(1) * AMBIENT_FRAC, soup)).astype(np.int32)
    E = r.poisson(np.outer(r.gamma(2.0, 20.0, N_EMPTY), soup)).astype(np.int32)

    genes = pd.DataFrame({"gene_symbol": [f"GENE{i}" for i in range(N_GENES)],
                          "feature_type": ["Gene Expression"]*N_GENES},
                         index=[f"ENSG{i:011d}" for i in range(N_GENES)])
    cell_bc  = [f"{name}_CELL{i:05d}" for i in range(X.shape[0])]
    empty_bc = [f"{name}_EMPTY{i:05d}" for i in range(N_EMPTY)]

    filt = ad.AnnData(sp.csr_matrix(X),
                      obs=pd.DataFrame({"sample": name, "condition": cond, "batch": name,
                                        "truth": labels}, index=cell_bc),
                      var=genes.copy())
    filt.layers["counts"] = filt.X.copy()
    raw = ad.AnnData(sp.csr_matrix(np.vstack([X, E])),
                     obs=pd.DataFrame({"sample": name, "condition": cond, "batch": name,
                                       "truth": labels + ["empty"]*N_EMPTY},
                                      index=cell_bc + empty_bc),
                     var=genes.copy())
    raw.layers["counts"] = raw.X.copy()
    filt.write_h5ad(f"h5ad/{name}.h5ad"); raw.write_h5ad(f"h5ad/{name}.raw.h5ad")
    return X.shape[0]

# FOUR samples, TWO per condition. With one sample per condition, batch and
# condition are the SAME variable and integrate.py correctly refuses to run:
# correcting the batch would delete the treatment effect. That is an
# experimental-design failure no setting can fix, so the fixture must encode a
# design that is actually analysable -- which is also the tutorial's design
# (4 healthy / 4 treated). Each sample still carries its own batch shift, so
# there is a real batch effect to remove WITHIN each condition.
specs = [("S1", 1.00, "healthy"), ("S2", 1.15, "healthy"),
         ("S3", 1.60, "treated"), ("S4", 1.75, "treated")]
counts = [make_sample(name, scale, 11 + i, cond) for i, (name, scale, cond) in enumerate(specs)]

pathlib.Path("ground_truth.json").write_text(json.dumps(
    {"n_samples": len(specs), "cells_per_sample": counts, "n_types": N_TYPES,
     "conditions": {c: [n for n, _, cc in specs if cc == c] for c in {x[2] for x in specs}},
     "doublets_per_sample": DOUBLET_N, "empties_per_sample": N_EMPTY,
     "n_genes": N_GENES, "n_marker_genes": N_MARKERS, "ambient_frac": AMBIENT_FRAC}, indent=2))
print(f"qc fixture: {len(specs)} samples x {counts} cells + {N_EMPTY} empties each, "
      f"{N_TYPES} types with EXCLUSIVE markers, {DOUBLET_N} doublets, "
      f"{AMBIENT_FRAC:.0%} ambient, batch effect planted")

# SELF-CHECK: does this fixture actually satisfy what SoupX needs?
chk = ad.read_h5ad("h5ad/S1.h5ad"); Xc = chk.X.toarray(); n = Xc.shape[0]
idf = np.log(n / np.maximum((Xc > 0).sum(0), 1))
m_idf, h_idf = idf[:MARKERS_PER_TYPE].mean(), idf[N_MARKERS:].mean()
print(f"  self-check: marker idf={m_idf:.2f}  housekeeping idf={h_idf:.2f}  "
      f"max idf={idf.max():.2f}")
if m_idf < 1.0:
    raise SystemExit(f"FIXTURE INVALID: marker idf {m_idf:.2f} < 1.0 -- ambient has made "
                     "markers detectable in nearly every cell, so tf-idf cannot "
                     "discriminate and SoupX will find no markers.\n"
                     "  idf <= log(N_TYPES) is a hard ceiling: a perfectly exclusive marker\n"
                     "  is still expressed in 1/N_TYPES of all cells. RAISE N_TYPES (currently "
                     f"{N_TYPES}, ceiling {np.log(N_TYPES):.2f}); lowering MARKER_LEVEL or\n"
                     "  AMBIENT_FRAC will barely move it.")
PY

cat > "$OUT/samplesheet.csv" <<'CSV'
sample,condition,batch
S1,healthy,S1
S2,healthy,S2
S3,treated,S3
S4,treated,S4
CSV
echo "qc fixture ready: $OUT"
