#!/usr/bin/env bash
# Install the R stack, conda-first with a BiocManager/CRAN source fallback.
# Same contract as install_quantifiers.sh: specs come from a yml so the pin list
# lives in one place, and every decision is echoed into the build log.
set -euo pipefail

ARCH="$(uname -m)"
SPEC_FILE="${1:-/tmp/environment.r.yml}"
PREFIX="${MAMBA_ROOT_PREFIX:-/opt/conda}"
echo "=== R stack install :: arch=${ARCH} :: specs from ${SPEC_FILE} ==="
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

# conda package name -> R package name, for the source fallback
r_name_for() {
  local n="$1"
  case "$n" in
    bioconductor-*) echo "${n#bioconductor-}" ;;
    r-base|r-biocmanager|r-remotes) echo "" ;;
    r-*) echo "${n#r-}" ;;
    *) echo "$n" ;;
  esac
}

mapfile -t SPECS < <(sed -n '/^dependencies:/,$p' "$SPEC_FILE" \
                     | grep -E '^[[:space:]]+-[[:space:]]+' \
                     | sed -E 's/^[[:space:]]*-[[:space:]]*//; s/[[:space:]]*#.*$//' \
                     | sed '/^$/d')

declare -a VIA_CONDA=() VIA_SOURCE=() FAILED=()

# r-base first: the source fallback needs a working R.
mm_install r-base=4.5 r-biocmanager r-remotes

# Installing package-by-package (below) lets every later transaction re-solve,
# and a Bioconductor package can silently BUMP r-base out from under us -- a
# pin declared in the yml but not enforced here does not hold. Pin it for real.
mkdir -p "${PREFIX:-/opt/conda}/conda-meta"
echo "r-base 4.5.*" >> "${PREFIX:-/opt/conda}/conda-meta/pinned"
echo ">>> pinned r-base 4.5.* for all subsequent transactions"

VIA_CONDA+=(r-base r-biocmanager r-remotes)

for spec in "${SPECS[@]}"; do
  name="${spec%%[<>=!]*}"
  case "$name" in r-base|r-biocmanager|r-remotes) continue ;; esac
  if conda_has "$spec"; then
    echo ">>> ${spec}: conda"
    mm_install "$spec"
    VIA_CONDA+=("$name")
  else
    rname="$(r_name_for "$name")"
    echo ">>> ${spec}: NO conda build for ${ARCH} -- BiocManager::install('${rname}')"
    if Rscript -e "BiocManager::install('${rname}', ask=FALSE, update=FALSE, Ncpus=$(nproc))" \
               -e "if (!requireNamespace('${rname}', quietly=TRUE)) quit(status=1)"; then
      VIA_SOURCE+=("$name")
    else
      echo "!!! ${name}: FAILED to install" >&2
      FAILED+=("$name")
    fi
  fi
done

# DoubletFinder is GitHub-only -- never on conda, so always this path.
Rscript -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', upgrade='never')" \
        -e "if (!requireNamespace('DoubletFinder', quietly=TRUE)) quit(status=1)" \
  && VIA_SOURCE+=(DoubletFinder) || FAILED+=(DoubletFinder)

micromamba clean -afy

echo "=== R stack install complete ==="
echo "    via conda : ${VIA_CONDA[*]:-none}"
echo "    via source: ${VIA_SOURCE[*]:-none}"
echo "    FAILED    : ${FAILED[*]:-none}"
echo "    r-base    : $(Rscript -e 'cat(paste(R.version$major, R.version$minor, sep="."))')"
[ ${#FAILED[@]} -eq 0 ] || { echo "R stack incomplete" >&2; exit 1; }
