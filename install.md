# Installation

The usual install is in [README.md](README.md): clone the repository, then `pixi install` and `pixi run setup-annot`. This file is for rebuilding the environment from scratch and for packages that pixi cannot always install from conda.

## Rebuild the pixi environment

Only needed if you are not using the committed `pixi.toml` / `pixi.lock` (for example when starting a new project).

```bash
## Intialize the pixi environment
pixi init
pixi workspace channel add bioconda
pixi project platform add osx-arm64 linux-64
```

Add the packages.

```bash
## Add packages
pixi add \
  r-optparse \
  r-readr \
  r-dplyr \
  r-tibble \
  r-tidyr \
  r-stringr \
  r-ggplot2 \
  r-ggridges \
  r-patchwork \
  r-gprofiler2 \
  r-complexupset \
  r-ggrepel \
  r-ggforce \
  r-tidydr \
  r-remotes \
  r-msigdbr \
  r-httpuv \
  r-shiny \
  r-xml2 \
  bioconductor-fgsea \
  bioconductor-complexheatmap \
  bioconductor-go.db \
  "bioconductor-org.hs.eg.db" \
  "bioconductor-org.mm.eg.db"
  ```

Activate the environment

```bash
pixi shell
```

## Annotation databases

`bioconductor-go.db`, `bioconductor-org.hs.eg.db`, and `bioconductor-org.mm.eg.db` are sub packages. Pixi does not run conda post-link scripts, so `pixi install` does not unpack the R databases. 

```bash
pixi run bash -c 'export PREFIX="$CONDA_PREFIX"
installBiocDataPackage.sh "go.db-3.22.0"
installBiocDataPackage.sh "org.hs.eg.db-3.22.0"
installBiocDataPackage.sh "org.mm.eg.db-3.22.0"'
```

## Installing Simona and simplifyEnrichment

```bash
pixi run Rscript -e 'install.packages(c("simona", "simplifyEnrichment"), repos = c("https://bioc.r-universe.dev", "https://cloud.r-project.org"))'
```

Note: `pixi add bioconductor-simona` is not used since that package has no `osx-arm64` conda build, and pixi must solve every platform listed in the project.

Running the above commands is equivalent to running the following from the repository root:

```bash
pixi run setup-annot
```
