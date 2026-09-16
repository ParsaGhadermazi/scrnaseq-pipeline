#!/usr/bin/env bash
# Build the synthetic 10x fixture used to verify the pipeline end to end.
#
# Ground truth: 40 cells x 5 genes, 20,000 UMIs.
# Barcodes are drawn from the REAL 10x v3 whitelist shipped inside the image,
# so bustools/alevin-fry barcode correction behaves as it would on real data.
#
#   usage: tests/make_fixture.sh <outdir> [image]
set -euo pipefail
OUT="${1:?usage: make_fixture.sh <outdir> [image]}"
IMG="${2:-parsaghadermazi/myrnaseqpipeline:0.2.0}"
# The image is local-only and never pushed, so a missing one makes docker try
# to pull from Docker Hub and fail with a confusing "access denied". Say what
# actually happened instead. (Docker prunes remove it -- it is unused between runs.)
if ! docker image inspect "$IMG" >/dev/null 2>&1; then
  echo "ERROR: image '$IMG' not found locally." >&2
  echo "  It is built, never pulled -- a 'docker image prune -a' removes it." >&2
  echo "  Rebuild with:" >&2
  echo "    docker buildx build --platform linux/arm64 -f containers/Dockerfile \\" >&2
  echo "      -t parsaghadermazi/myrnaseqpipeline:0.2.0 --load ." >&2
  exit 1
fi

mkdir -p "$OUT"/{data,index}
OUT="$(cd "$OUT" && pwd)"

docker run --rm -i -v "$OUT":/w -w /w "$IMG" python - <<'PY'
import gzip, random, pathlib
random.seed(7)
WL = "/opt/conda/lib/python3.11/site-packages/ngs_tools/chemistry/whitelists/10x_version3_whitelist.txt.gz"
with gzip.open(WL, "rt") as f:
    cells = random.sample([next(f).strip() for _ in range(200)], 40)

tx, t2g = [], []
for g in range(5):
    for t in range(4):
        tid, gid = f"ENST{g:05d}{t}", f"ENSG{g:011d}"
        tx.append((tid, "".join(random.choice("ACGT") for _ in range(1500))))
        t2g.append(f"{tid}\t{gid}\tGENE{g}")
pathlib.Path("index/tx.fa").write_text("".join(f">{i}\n{s}\n" for i, s in tx))
pathlib.Path("index/t2g.txt").write_text("\n".join(t2g) + "\n")

r1, r2 = [], []
for i in range(20000):
    bc, umi = random.choice(cells), "".join(random.choice("ACGT") for _ in range(12))
    _, seq = random.choice(tx); st = random.randint(0, len(seq) - 91)
    r1.append(f"@r{i}\n{bc}{umi}\n+\n{'I'*28}\n")          # 28bp = 16 barcode + 12 UMI
    r2.append(f"@r{i}\n{seq[st:st+91]}\n+\n{'I'*91}\n")    # 91bp cDNA
for fn, rec in (("data/SRRTEST_1.fastq.gz", r1), ("data/SRRTEST_2.fastq.gz", r2)):
    with gzip.open(fn, "wt") as f: f.write("".join(rec))
print("fixture: 40 cells x 5 genes, 20000 read pairs")
PY

docker run --rm -v "$OUT":/w -w /w "$IMG" kallisto index -i index/index.idx index/tx.fa >/dev/null 2>&1
docker run --rm -e ALEVIN_FRY_HOME=/w/afhome -v "$OUT":/w -w /w "$IMG" \
  bash -lc 'simpleaf set-paths >/dev/null 2>&1; simpleaf index --output idx_af --ref-seq index/tx.fa -t 4' >/dev/null 2>&1
awk '{print $1"\t"$2}' "$OUT/index/t2g.txt" > "$OUT/idx_af/t2g.tsv"

printf 'sample,fastq_1,fastq_2,condition,batch\nSRRTEST,%s/data/SRRTEST_1.fastq.gz,%s/data/SRRTEST_2.fastq.gz,healthy,run1\n' \
  "$OUT" "$OUT" > "$OUT/samplesheet.csv"
echo "fixture ready: $OUT"
