# enrichviz

`enrichviz` provides R scripts for the visualization of functional enrichment analysis results.

It supports two complementary approaches:

- **Over-representation analysis (ORA):** identifies functions that occur more often among significant genes than expected by chance.
- **Gene set enrichment analysis (GSEA):** identifies functions whose genes are concentrated toward the top or bottom of a ranked list of all tested genes.

Both ORA and GSEA scripts support **GO** and **KEGG** results and generate publication-ready PDF plots together with the tables used to create them.

Detailed guides:

- [ORA visualization guide](ora/README_ORA.md)
- [GSEA visualization guide](gsea/README_GSEA.md)

## Installation

The environment is managed with [pixi](https://pixi.sh/) using `pixi.toml` and `pixi.lock`.

```bash
git clone https://github.com/aswathyseb/enrichviz.git
cd enrichviz
pixi install
pixi run setup-annot
```

`setup-annot` installs `go.db`, the human and mouse OrgDb packages, and `simona` / `simplifyEnrichment` (needed for `ora simplify`). The annotation download is required because `pixi install` does not run conda post-link scripts. 

Rebuilding the environment from scratch is documented in [install.md](install.md).

## Running enrichviz

The `enrichviz` command provides a common interface to the ORA and GSEA scripts. 

Start a pixi shell:

```bash
pixi shell
```

Then run commands directly. To get help run

```bash
enrichviz --help
```

Some example commands are below.

Perform ORA on GO terms

```bash
enrichviz ora gprofiler --in edger.csv --outdir res --direction yes
```
Create barplot of top N terms

```bash
enrichviz ora barplot --in res/gprofiler_GO.csv --outdir res --ont BP
```

Get help on individual commands with `-h`.

```bash
enrichviz ora barplot -h
```

Tab completion is available inside `pixi shell` for tracks, subcommands, and file paths.


## Available commands

| Command | Purpose |
|---|---|
| `enrichviz ora gprofiler` | Run GO enrichment with gProfiler |
| `enrichviz ora barplot` | ORA barplot |
| `enrichviz ora lollipop` | ORA lollipop plot |
| `enrichviz ora dotplot` | ORA dotplot |
| `enrichviz ora upset` | UpSet plot |
| `enrichviz ora ssplot` | Semantic space plot |
| `enrichviz ora simplify` | Simplify GO terms |
| `enrichviz gsea go` | Run GO GSEA |
| `enrichviz gsea kegg` | Run KEGG GSEA |
| `enrichviz gsea barplot` | GSEA barplot |
| `enrichviz gsea lollipop` | GSEA lollipop plot |
| `enrichviz gsea dotplot` | GSEA dotplot |
| `enrichviz gsea ridgeplot` | GSEA ridgeplot |
| `enrichviz gsea volcano` | GSEA volcano plot |
| `enrichviz gsea nes` | Running ES/NES plot |
| `enrichviz gsea nes-compare` | Compare NES across GSEA runs |

## ORA or GSEA?

ORA and GSEA answer related but different biological questions.

| | ORA | GSEA |
|---|---|---|
| **Input** | Significant DEGs and a background gene set | All tested genes ranked by a statistic |
| **Question** | Which functions are over-represented among DEGs? | Which functions are enriched toward either end of the ranked gene list? |
| **Direction** | Up- and down-regulated genes can be tested separately | Direction is represented by the ranking statistic and NES |
| **Run** | `enrichviz ora gprofiler` | `enrichviz gsea go` or `enrichviz gsea kegg` |

Use **ORA** when you have a defined list of significant DEGs and want to identify functions enriched within that list.

Use **GSEA** when you want to analyze the complete ranked gene list without applying a significance cutoff.

Plotting commands accept CSV, TSV, or RDS results from **gProfiler**, **clusterProfiler**, or **fgsea**.

## Choosing a visualization

Different plots highlight different aspects of the enrichment results.

**Summarizing the strongest enriched terms**

- **Barplot:** shows the magnitude and direction of the strongest terms.
- **Lollipop plot:** provides a similar overview while using point size to represent gene-set size or gene count.
- **Dotplot:** shows the fraction of genes contributing to enrichment.

For ORA, the dotplot represents the query gene ratio. For GSEA, it represents the leading-edge ratio.

**Exploring relationships among ORA terms**

- **UpSet plot:** shows genes shared among selected enriched terms.
- **Semantic space plot:** arranges terms according to the overlap of their query genes.
- **simplifyGO:** groups GO terms according to semantic similarity within the ontology.

**Examining GSEA signals**

- **Ridgeplot:** shows the distribution of ranking scores for leading-edge genes.
- **Volcano plot:** shows NES against adjusted *P* value for all tested terms.
- **Running ES/NES plot:** shows the enrichment trajectory of an individual term.
- **NES comparison:** compares the same term across multiple GSEA analyses.

For plot-specific options, term-selection rules, and example commands, see:

- [ORA visualization guide](ora/README_ORA.md)
- [GSEA visualization guide](gsea/README_GSEA.md)

## License

Released under the [MIT License](LICENSE).