#!/usr/bin/env bash
# End-to-end verification of stage `qc` against the synthetic fixture.
#   usage: tests/run_qc_tests.sh [workdir] [image]
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:-/tmp/scrnaseq_qc_tests}"
IMG="${2:-parsaghadermazi/myrnaseqpipeline:0.2.0}"
PASS=0; FAIL=0
ok(){ echo "  PASS  $1"; PASS=$((PASS+1)); }
no(){ echo "  FAIL  $1"; FAIL=$((FAIL+1)); }

echo "== fixture =="
bash "$HERE/tests/make_qc_fixture.sh" "$WORK" "$IMG" || { echo "fixture failed"; exit 1; }

echo "== running stage qc =="
# -profile docker,test : the `test` profile carries soupx_soup_quantile=0.50,
# without which SoupX fails the soup gate on a synthetic fixture (see docs/testing.md).
( cd "$WORK" && nextflow run "$HERE/main.nf" -profile docker,test \
    --stage qc --input "$WORK/samplesheet.csv" --outdir "$WORK/out" \
    --counts_dir "$WORK/h5ad" --max_cpus 4 --max_memory 6.GB ) 2>&1 | tail -20

I="$WORK/out/qc/integrated/integrated.h5ad"
[ -f "$I" ] && ok "integrated.h5ad produced" || { no "integrated.h5ad missing"; echo "== $PASS passed, $FAIL failed =="; exit 1; }

echo "== output contract + did integration actually do anything? =="
cat > "$WORK/_qccheck.py" <<'PYCHECK'
import anndata as ad, numpy as np, json, sys
a = ad.read_h5ad('out/qc/integrated/integrated.h5ad')
t = json.load(open('ground_truth.json'))
planted = sum(t['cells_per_sample'])

c = [
    ('counts layer preserved (DE needs raw)', 'counts' in a.layers),
    ('X_integrated present',                  'X_integrated' in a.obsm),
    ('X_pca retained (uncorrected)',          'X_pca' in a.obsm),
    ('both samples present',                  a.obs['sample'].nunique() == t['n_samples']),
    ('cells reduced by QC',                   a.n_obs < planted),
    ('cells not wiped out',                   a.n_obs > 0.4 * planted),
    ('integration recorded in uns',           'integration' in a.uns),
]

# Structural checks alone cannot tell real integration from a no-op copy of
# X_pca. These two can.
if 'X_integrated' in a.obsm and 'X_pca' in a.obsm:
    P, I = np.asarray(a.obsm['X_pca']), np.asarray(a.obsm['X_integrated'])
    same_shape = P.shape == I.shape
    moved = (not same_shape) or not np.allclose(P, I)
    c.append(('X_integrated is NOT a copy of X_pca', moved))

    # the fixture plants a batch effect, so correction must REDUCE separability
    # of samples in the embedding -- measured as batch silhouette
    from sklearn.metrics import silhouette_score
    b = a.obs['sample'].astype(str).to_numpy()
    if len(np.unique(b)) > 1 and same_shape:
        n = min(2000, a.n_obs)
        idx = np.random.default_rng(0).choice(a.n_obs, n, replace=False)
        s_pre  = silhouette_score(P[idx], b[idx])
        s_post = silhouette_score(I[idx], b[idx])
        print(f"  batch silhouette: uncorrected {s_pre:+.4f} -> integrated {s_post:+.4f}")
        c.append(('batch silhouette DECREASED (batch effect removed)', s_post < s_pre))

for name, good in c:
    print(('  PASS  ' if good else '  FAIL  ') + name)
sys.exit(0 if all(g for _, g in c) else 1)
PYCHECK
docker run --rm -v "$WORK":/w -w /w "$IMG" python /w/_qccheck.py 2>&1 | grep -vi "futurewarning"
rc=${PIPESTATUS[0]}
[ $rc -eq 0 ] && PASS=$((PASS+9)) || FAIL=$((FAIL+1))

M="$WORK/out/qc/integrated/integration_metrics.json"
if [ -f "$M" ]; then
  python3 -c "
import json; d=json.load(open('$M')); m=d.get('metrics',{})
print('  metrics:', json.dumps(d.get('delta', m)))
assert 'integrated' in m, 'no integrated metrics'
" && ok "integration metrics emitted" || no "integration metrics malformed"
else no "integration_metrics.json missing"; fi

echo; echo "== $PASS passed, $FAIL failed =="
exit $([ $FAIL -eq 0 ] && echo 0 || echo 1)
