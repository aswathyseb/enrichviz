# Installation 

## Initialize pixi environment
pixi init

## Add the bioconda channel
pixi workspace channel add bioconda

## The platforms that the environment will support
pixi project platform add osx-arm64 linux-64

## Activate the environment

Activate the environment with the command:

```bash
pixi shell
```

## Add packages

```bash
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
  bioconductor-fgsea \
  bioconductor-complexheatmap \
  bioconductor-go.db \
  "bioconductor-org.hs.eg.db" \
  "bioconductor-org.mm.eg.db"
```

## Download annotation databases

`bioconductor-go.db`, `bioconductor-org.hs.eg.db`, and `bioconductor-org.mm.eg.db` are stub packages. Pixi does not run conda post-link scripts, so `pixi add` does not install the R packages. Download them into the pixi library:

```bash
pixi run bash -c 'export PREFIX="$CONDA_PREFIX"
installBiocDataPackage.sh "go.db-3.22.0"
installBiocDataPackage.sh "org.hs.eg.db-3.22.0"
installBiocDataPackage.sh "org.mm.eg.db-3.22.0"'
```

## Install simona and simplifyEnrichment

`bioconductor-simona` has no `osx-arm64` conda build, so `pixi add bioconductor-simona` cannot be solved on Apple Silicon. `simplifyEnrichment` depends on `simona`.

On linux-64 or Intel macOS:

```bash
pixi add bioconductor-simona bioconductor-simplifyenrichment
```

On Apple Silicon (`osx-arm64`): 

```bash
pixi run Rscript -e 'install.packages(c("simona", "simplifyEnrichment"), repos = c("https://bioc-release.r-universe.dev", "https://cloud.r-project.org"))'
```






