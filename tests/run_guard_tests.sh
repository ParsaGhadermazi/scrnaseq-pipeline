#!/usr/bin/env bash
# Negative tests: the pipeline must REFUSE bad input rather than produce
# plausible-looking wrong output.
#
# These guard the failure mode that actually matters. A crash is obvious; a run
# that silently removes the treatment effect and emits a clean-looking object is
# the one that reaches a paper.
#
#   usage: tests/run_guard_tests.sh [workdir] [image]
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:-/tmp/scrnaseq_guard_tests}"
IMG="${2:-parsaghadermazi/myrnaseqpipeline:0.2.0}"
PASS=0; FAIL=0
ok(){ echo "  PASS  $1"; PASS=$((PASS+1)); }
no(){ echo "  FAIL  $1"; FAIL=$((FAIL+1)); }

echo "== fixture (reused for all guard cases) =="
[ -f "$WORK/h5ad/S1.h5ad" ] || bash "$HERE/tests/make_qc_fixture.sh" "$WORK" "$IMG" >/dev/null \
  || { echo "fixture failed"; exit 1; }

# --- 1. batch perfectly confounded with condition -------------------------
# One sample per condition => batch IS condition. Correcting the batch would
# delete the treatment effect, and no parameter can fix it.
cat > "$WORK/ss_confounded.csv" <<'CSV'
sample,condition,batch
S1,healthy,S1
S3,treated,S3
CSV
out=$( cd "$WORK" && nextflow run "$HERE/main.nf" -profile docker,test \
        --stage qc --input "$WORK/ss_confounded.csv" --outdir "$WORK/out_conf" \
        --counts_dir "$WORK/h5ad" --max_cpus 2 --max_memory 4.GB 2>&1 )
echo "$out" | grep -q "perfectly confounded" \
  && ok "confounded batch/condition refused with an explanatory error" \
  || no "confounded design NOT refused (would silently delete the condition effect)"

# --- 2. missing raw matrix -------------------------------------------------
# SoupX and emptyDrops measure the soup FROM the empty droplets, which exist
# only in the raw matrix. Starting without it must fail loudly, not silently
# skip ambient correction.
mkdir -p "$WORK/h5ad_noraw" && cp "$WORK/h5ad/S1.h5ad" "$WORK/h5ad_noraw/" 2>/dev/null
printf 'sample,condition,batch\nS1,healthy,S1\n' > "$WORK/ss_noraw.csv"
out=$( cd "$WORK" && nextflow run "$HERE/main.nf" -profile docker,test \
        --stage qc --input "$WORK/ss_noraw.csv" --outdir "$WORK/out_noraw" \
        --counts_dir "$WORK/h5ad_noraw" --max_cpus 2 --max_memory 4.GB 2>&1 )
echo "$out" | grep -qi "raw" \
  && ok "missing .raw.h5ad refused" \
  || no "missing raw matrix NOT refused"

echo; echo "== $PASS passed, $FAIL failed =="
exit $([ $FAIL -eq 0 ] && echo 0 || echo 1)
