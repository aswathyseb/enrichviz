#!/usr/bin/env bash
# Download Bioconductor annotation databases. Pixi/conda stub packages do not
# run post-link scripts, so go.db and OrgDb must be fetched explicitly.
set -euo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "CONDA_PREFIX is not set. Run: pixi run setup-annot" >&2
  exit 1
fi

if ! command -v installBiocDataPackage.sh >/dev/null 2>&1; then
  echo "installBiocDataPackage.sh not found. Run pixi install first." >&2
  exit 1
fi

export PREFIX="$CONDA_PREFIX"
installBiocDataPackage.sh "go.db-3.22.0"
installBiocDataPackage.sh "org.hs.eg.db-3.22.0"
installBiocDataPackage.sh "org.mm.eg.db-3.22.0"
