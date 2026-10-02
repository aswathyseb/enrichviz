# Download GO.db and the human/mouse OrgDb packages, then install simona and
# simplifyEnrichment. Run with: pixi run setup-annot
#
# The conda builds of these annotation packages are stubs. This script
# downloads the Bioconductor 3.22 tarballs and installs them into the pixi
# library. simona has no osx-arm64 conda build, so it and simplifyEnrichment
# come from r-universe.

prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
if (!nzchar(prefix)) {
  stop("CONDA_PREFIX is not set. Run: pixi run setup-annot", call. = FALSE)
}

lib <- .libPaths()[[1]]
if (!startsWith(normalizePath(lib, winslash = "/"), normalizePath(prefix, winslash = "/"))) {
  stop("R library is not inside CONDA_PREFIX: ", lib, call. = FALSE)
}

options(timeout = max(600, getOption("timeout")))

annot_packages <- list(
  list(
    pkg = "GO.db",
    sqlite = "GO.sqlite",
    md5 = "5ae5557afa56227c4c9c145907b1f585",
    urls = c(
      "https://bioconductor.org/packages/3.22/data/annotation/src/contrib/GO.db_3.22.0.tar.gz",
      "https://bioarchive.galaxyproject.org/GO.db_3.22.0.tar.gz",
      "https://depot.galaxyproject.org/software/bioconductor-go.db/bioconductor-go.db_3.22.0_src_all.tar.gz"
    )
  ),
  list(
    pkg = "org.Hs.eg.db",
    sqlite = "org.Hs.eg.sqlite",
    md5 = "e80cac6ec018a95aea4f7530350e80a2",
    urls = c(
      "https://bioconductor.org/packages/3.22/data/annotation/src/contrib/org.Hs.eg.db_3.22.0.tar.gz",
      "https://bioarchive.galaxyproject.org/org.Hs.eg.db_3.22.0.tar.gz",
      "https://depot.galaxyproject.org/software/bioconductor-org.hs.eg.db/bioconductor-org.hs.eg.db_3.22.0_src_all.tar.gz"
    )
  ),
  list(
    pkg = "org.Mm.eg.db",
    sqlite = "org.Mm.eg.sqlite",
    md5 = "cc19f63b25a7a8018e46f71e42e2f83f",
    urls = c(
      "https://bioconductor.org/packages/3.22/data/annotation/src/contrib/org.Mm.eg.db_3.22.0.tar.gz",
      "https://bioarchive.galaxyproject.org/org.Mm.eg.db_3.22.0.tar.gz",
      "https://depot.galaxyproject.org/software/bioconductor-org.mm.eg.db/bioconductor-org.mm.eg.db_3.22.0_src_all.tar.gz"
    )
  )
)

install_annot <- function(spec, lib) {
  db <- file.path(lib, spec$pkg, "extdata", spec$sqlite)
  if (file.exists(db)) {
    message(spec$pkg, " already installed")
    return(invisible())
  }

  dest <- file.path(tempdir(), paste0(spec$pkg, "_3.22.0.tar.gz"))
  ok <- FALSE
  for (url in spec$urls) {
    message("Downloading ", url)
    status <- tryCatch(
      download.file(url, dest, mode = "wb", quiet = FALSE),
      error = function(e) 1L
    )
    if (!is.numeric(status) || status != 0L || !file.exists(dest)) {
      next
    }
    got <- unname(tools::md5sum(dest))
    if (identical(got, spec$md5)) {
      ok <- TRUE
      break
    }
    message("Checksum mismatch for ", url)
    unlink(dest)
  }
  if (!ok) {
    stop("Failed to download ", spec$pkg, call. = FALSE)
  }

  install.packages(dest, repos = NULL, type = "source", lib = lib)
  unlink(dest)
  if (!file.exists(db)) {
    stop("Installed ", spec$pkg, " but ", spec$sqlite, " is missing", call. = FALSE)
  }
  invisible()
}

for (spec in annot_packages) {
  install_annot(spec, lib)
}

need <- c("simona", "simplifyEnrichment")
need <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(need) == 0L) {
  message("simona and simplifyEnrichment already installed")
} else {
  install.packages(
    need,
    lib = lib,
    repos = c(
      "https://bioc.r-universe.dev",
      "https://cloud.r-project.org"
    )
  )
  missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Failed to install: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}
