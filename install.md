# Installation

The usual install is in [README.md](README.md): clone the repository, then `pixi install` and `pixi run setup-annot`. This file is for rebuilding the environment from scratch and for packages that pixi cannot always install from conda.

## Rebuild the pixi environment

Only needed if you are not using the committed `pixi.toml` / `pixi.lock` (for example when starting a new project).

```bash
pixi init
pixi workspace channel add bioconda
pixi project platform add osx-arm64 linux-64
```

Then add the packages listed in `pixi.toml`, or copy that file and run `pixi install`.

Activate the environment when you want an interactive prompt (`R`, `enrichviz` on `PATH`):

```bash
pixi shell
```

## Annotation databases

`bioconductor-go.db`, `bioconductor-org.hs.eg.db`, and `bioconductor-org.mm.eg.db` are stub packages. Pixi does not run conda post-link scripts, so `pixi install` does not unpack the R databases. From the repository root:

```bash
pixi run setup-annot
```

That runs `scripts/install-annot.sh`, which is equivalent to:

```bash
pixi run bash -c 'export PREFIX="$CONDA_PREFIX"
installBiocDataPackage.sh "go.db-3.22.0"
installBiocDataPackage.sh "org.hs.eg.db-3.22.0"
installBiocDataPackage.sh "org.mm.eg.db-3.22.0"'
```

## Install simona and simplifyEnrichment

`enrichviz ora simplify` needs these packages. `bioconductor-simona` has no `osx-arm64` conda build, so `pixi add bioconductor-simona` cannot be solved on Apple Silicon. `simplifyEnrichment` depends on `simona`.

On linux-64 or Intel macOS:

```bash
pixi add bioconductor-simona bioconductor-simplifyenrichment
```

On Apple Silicon (`osx-arm64`):

```bash
pixi run Rscript -e 'install.packages(c("simona", "simplifyEnrichment"), repos = c("https://bioc-release.r-universe.dev", "https://cloud.r-project.org"))'
```
