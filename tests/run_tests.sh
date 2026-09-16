#!/usr/bin/env bash
# End-to-end verification against synthetic data with known ground truth.
#
#   usage: tests/run_tests.sh [workdir]
#
# Ground truth: 40 cells x 5 genes, 20,000 UMIs.
# Checks equality against those numbers rather than eyeballing output.
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:-/tmp/scrnaseq_tests}"
IMG=parsaghadermazi/myrnaseqpipeline:0.2.0
PASS=0; FAIL=0
ok(){ echo "  PASS  $1"; PASS=$((PASS+1)); }
no(){ echo "  FAIL  $1"; FAIL=$((FAIL+1)); }

echo "== building fixture =="
bash "$HERE/tests/make_fixture.sh" "$WORK" >/dev/null || { echo "fixture build failed"; exit 1; }

echo "== 1. detector: chemistry mis-declaration must be refused =="
docker run --rm -v "$WORK":/w -v "$HERE/bin":/b -w /w $IMG \
  python /b/detect_10x_reads.py --fastqs data/SRRTEST_1.fastq.gz data/SRRTEST_2.fastq.gz \
  --sample T --outdir /tmp/x --json /tmp/x.json --declared-chemistry SC3Pv2 >/dev/null 2>&1
[ $? -ne 0 ] && ok "mis-declared chemistry rejected" || no "mis-declared chemistry accepted"

echo "== 2. detector: missing barcode read must be refused =="
docker run --rm -v "$WORK":/w -v "$HERE/bin":/b -w /w $IMG \
  python /b/detect_10x_reads.py --fastqs data/SRRTEST_2.fastq.gz \
  --sample T --outdir /tmp/y --json /tmp/y.json >/dev/null 2>&1
[ $? -ne 0 ] && ok "missing barcode read rejected" || no "missing barcode read accepted"

echo "== 3. both quantifiers, end to end =="
for spec in "kallisto:$WORK/index" "alevinfry:$WORK/idx_af"; do
  q=${spec%%:*}; idx=${spec##*:}
  ( cd "$WORK" && nextflow run "$HERE/main.nf" -profile docker \
      --input "$WORK/samplesheet.csv" --outdir "$WORK/out_$q" \
      --quantifier "$q" --index "$idx" --max_cpus 4 --max_memory 6.GB ) >/dev/null 2>&1
  [ -f "$WORK/out_$q/h5ad/SRRTEST.h5ad" ] && ok "$q produced h5ad" || no "$q produced no h5ad"
done

echo "== 4. ground truth + cross-quantifier schema parity =="
cat > "$WORK/_check.py" <<'PYCHECK'
import anndata as ad, numpy as np, sys
k = ad.read_h5ad('out_kallisto/h5ad/SRRTEST.h5ad')
a = ad.read_h5ad('out_alevinfry/h5ad/SRRTEST.h5ad')
checks = [
    ('kallisto 40x5',      k.shape == (40, 5)),
    ('alevinfry 40x5',     a.shape == (40, 5)),
    ('kallisto sum 20000', int(k.X.sum()) == 20000),
    ('obs schema equal',   list(k.obs.columns) == list(a.obs.columns)),
    ('var schema equal',   list(k.var.columns) == list(a.var.columns)),
    ('gene sets equal',    set(k.var_names) == set(a.var_names)),
]
g = sorted(set(k.var_names) & set(a.var_names))
o = sorted(set(k.obs_names) & set(a.obs_names))
r = np.corrcoef(k[o, g].X.toarray().ravel(), a[o, g].X.toarray().ravel())[0, 1]
checks.append((f'concordance r={r:.4f} > 0.99', r > 0.99))
for name, good in checks:
    print(('  PASS  ' if good else '  FAIL  ') + name)
sys.exit(0 if all(good for _, good in checks) else 1)
PYCHECK
docker run --rm -v "$WORK":/w -w /w "$IMG" python /w/_check.py 2>&1 | grep -vi "futurewarning"
rc=${PIPESTATUS[0]}
[ $rc -eq 0 ] && PASS=$((PASS+7)) || FAIL=$((FAIL+1))

echo
echo "== $PASS passed, $FAIL failed =="
exit $([ $FAIL -eq 0 ] && echo 0 || echo 1)
