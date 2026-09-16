#!/usr/bin/env bash
# Download Bioconductor annotation databases and install simona /
# simplifyEnrichment. Pixi/conda stub packages do not run post-link scripts,
# so go.db and OrgDb must be fetched explicitly. simona has no osx-arm64
# conda build, so those two R packages come from r-universe.
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

Rscript -e '
need <- c("simona", "simplifyEnrichment")
need <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(need) == 0L) {
  message("simona and simplifyEnrichment already installed")
  quit(save = "no", status = 0)
}
install.packages(
  need,
  repos = c("https://bioc.r-universe.dev", "https://cloud.r-project.org")
)
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Failed to install: ", paste(missing, collapse = ", "), call. = FALSE)
}
'
