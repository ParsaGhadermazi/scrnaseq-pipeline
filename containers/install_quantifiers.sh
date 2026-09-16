#!/usr/bin/env bash
# Install the quantifier stack, preferring conda and falling back to a source
# build only when the package genuinely has no build for the target arch.
#
# Package specs are read from env/environment.quant.yml so the pin list lives in
# exactly one place. Every decision is echoed, so the build log is the record of
# what actually happened on this architecture.
set -euo pipefail

ARCH="$(uname -m)"
PREFIX="${MAMBA_ROOT_PREFIX:-/opt/conda}"
SPEC_FILE="${1:-/tmp/environment.quant.yml}"
SRC=/tmp/src
mkdir -p "$SRC"
echo "=== quantifier install :: arch=${ARCH} :: specs from ${SPEC_FILE} ==="

CHANNELS=(-c conda-forge -c bioconda)

# bioconda repodata fetches time out intermittently ("Download error (28)").
# That is transient network, not a solver problem, so retry rather than failing
# a 20-minute build over it.
mm_install() {
  local n=0
  until micromamba install -y -n base "${CHANNELS[@]}" "$@"; do
    n=$((n+1))
    [ $n -ge 3 ] && { echo "!!! micromamba install failed 3x: $*" >&2; return 1; }
    echo "    retry $n/3 after repodata/network failure: $*" >&2
    sleep $((n*20))
  done
}

# NOTE: the channels are essential here. Without them this resolves against the
# default channels only and reports false negatives for every bioconda package,
# sending everything down the source-build path.
# A dry-run that FAILS is ambiguous: the package may genuinely not exist for
# this arch, or bioconda's repodata may just have timed out. Treating both as
# "not available" silently reroutes a package to a source build and changes what
# ends up in the image -- which is exactly what happened to sra-tools on a
# flaky-network build. So classify the failure and abort on network errors
# rather than quietly building from source.
conda_has() {
  local spec="$1" out rc n=0
  while :; do
    out="$(micromamba install -y -n base "${CHANNELS[@]}" --dry-run "$spec" 2>&1)"; rc=$?
    [ $rc -eq 0 ] && return 0
    if echo "$out" | grep -qiE 'download error|timeout|too slow|repodata|could not resolve|connection'; then
      n=$((n+1))
      if [ $n -ge 3 ]; then
        echo "!!! repodata unreachable while probing '${spec}' after 3 tries." >&2
        echo "    Refusing to guess availability -- a network failure must not" >&2
        echo "    silently become a source build. Re-run when the network is up." >&2
        exit 1
      fi
      echo "    repodata probe failed for '${spec}', retry ${n}/3" >&2
      sleep $((n*20)); continue
    fi
    return 1          # genuine "no build for this arch"
  done
}

# --- source-build fallbacks -------------------------------------------------

build_rust() {   # Rust; builds cleanly on both arches
  command -v cargo >/dev/null || { curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal; }
  . "$HOME/.cargo/env" 2>/dev/null || true
  cargo install "$1" --root "$PREFIX" --locked
}

build_cmake() {  # $1=name $2=url $3=tag
  git clone --depth 1 --branch "$3" --recursive "$2" "$SRC/$1"
  cmake -S "$SRC/$1" -B "$SRC/$1/build" -DCMAKE_INSTALL_PREFIX="$PREFIX" -DCMAKE_BUILD_TYPE=Release
  cmake --build "$SRC/$1/build" -j"$(nproc)" --target install
}

build_sratools() {
  if [ "$ARCH" = "x86_64" ]; then
    curl -sSL https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz \
      | tar xz -C "$SRC"
    cp "$SRC"/sratoolkit.*/bin/* "$PREFIX/bin/" 2>/dev/null || true
    return
  fi
  git clone --depth 1 https://github.com/ncbi/ncbi-vdb "$SRC/ncbi-vdb"
  cmake -S "$SRC/ncbi-vdb" -B "$SRC/ncbi-vdb/build" -DCMAKE_INSTALL_PREFIX="$PREFIX"
  cmake --build "$SRC/ncbi-vdb/build" -j"$(nproc)" --target install
  # ncbi-vdb installs into lib64/ but sra-tools' cmake looks in lib/.
  # Bridge them rather than guessing which one a given version picks.
  for f in "$PREFIX"/lib64/libncbi-vdb*; do
    [ -e "$f" ] && ln -sf "$f" "$PREFIX/lib/$(basename "$f")"
  done
  git clone --depth 1 https://github.com/ncbi/sra-tools "$SRC/sra-tools"
  cmake -S "$SRC/sra-tools" -B "$SRC/sra-tools/build" -DCMAKE_INSTALL_PREFIX="$PREFIX" \
        -DVDB_LIBDIR="$PREFIX/lib" -DVDB_INCDIR="$PREFIX/include"
  cmake --build "$SRC/sra-tools/build" -j"$(nproc)" --target install
}

build_pip() { pip install --no-cache-dir "$1"; }

fallback_for() {  # package name -> source-build command
  case "$1" in
    sra-tools)  build_sratools ;;
    salmon)     build_cmake salmon   https://github.com/COMBINE-lab/salmon v1.10.3 ;;
    kallisto)   build_cmake kallisto https://github.com/pachterlab/kallisto v0.50.1 ;;
    bustools)   build_cmake bustools https://github.com/BUStools/bustools v0.43.2 ;;
    alevin-fry) build_rust alevin-fry ;;
    simpleaf)   build_rust simpleaf ;;
    *)          build_pip "$2" ;;      # kb-python, pyroe, anything new
  esac
}

# --- drive from the spec file ------------------------------------------------

mapfile -t SPECS < <(sed -n '/^dependencies:/,$p' "$SPEC_FILE" \
                     | grep -E '^[[:space:]]+-[[:space:]]+' \
                     | sed -E 's/^[[:space:]]*-[[:space:]]*//; s/[[:space:]]*#.*$//' \
                     | sed '/^$/d')

echo "specs: ${SPECS[*]}"
declare -a VIA_CONDA=() VIA_SOURCE=()

for spec in "${SPECS[@]}"; do
  name="${spec%%[<>=!]*}"
  if conda_has "$spec"; then
    echo ">>> ${spec}: conda"
    mm_install "$spec"
    VIA_CONDA+=("$name")
  else
    echo ">>> ${spec}: NO conda build for ${ARCH} -- building from source"
    fallback_for "$name" "$spec"
    VIA_SOURCE+=("$name")
  fi
done

micromamba clean -afy
rm -rf "$SRC" "$HOME/.cargo/registry"

echo "=== quantifier install complete ==="
echo "    via conda : ${VIA_CONDA[*]:-none}"
echo "    via source: ${VIA_SOURCE[*]:-none}"
