# Installation

The usual install is explained in [README.md](README.md) where the steps are
1. Clone the repository
2. Run the command `pixi install` 
3. Run  `pixi run setup-annot`. 

Described below are the steps to built the environment from scratch. These are needed only if you are not using the provided `pixi.toml` and `pixi.lock` 

## Rebuild the pixi environment

Only needed if you are not using the committed `pixi.toml` / `pixi.lock`.

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
  bioconductor-clusterprofiler \
  bioconductor-annotationdbi \
  bioconductor-go.db \
  "bioconductor-org.hs.eg.db" \
  "bioconductor-org.mm.eg.db"
  ```

Activate the environment

```bash
pixi shell
```

## Annotation packages

`pixi install` does not unpack the conda stub packages `bioconductor-go.db`, `bioconductor-org.hs.eg.db`, and `bioconductor-org.mm.eg.db`. `pixi run setup-annot` runs `src/install-annot.R`, which downloads the Bioconductor 3.22 tarballs for `GO.db`, `org.Hs.eg.db`, and `org.Mm.eg.db` and installs them into the pixi library.

The same script installs `simona` and `simplifyEnrichment` from r-universe. `pixi add bioconductor-simona` is not used, because that package has no `osx-arm64` conda build and pixi must solve every platform listed in the project.

Running all the above commands is equivalent to running the following from the repository root:

```bash
pixi install
pixi run setup-annot
```
