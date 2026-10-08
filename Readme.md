# enrichviz

`enrichviz` provides R scripts for the visualization of functional enrichment analysis results.

It supports two complementary approaches and generate publication-ready PDF plots.


- **Over-representation analysis (ORA):** identifies functions that occur more often among significant genes than expected by chance.
- **Gene set enrichment analysis (GSEA):** identifies functions whose genes are concentrated toward the top or bottom of a ranked list of all tested genes.

Both ORA and GSEA scripts support **GO** and **KEGG** based 

Detailed guides:

- [ORA visualization guide](docs/ora.md)
- [GSEA visualization guide](docs/gsea.md)

## Installation

The environment is managed with [pixi](https://pixi.sh/) using `pixi.toml` and `pixi.lock`.

```bash
git clone https://github.com/aswathyseb/enrichviz.git
cd enrichviz
pixi install
pixi run setup-annot
```

**Note:**
`setup-annot` installs  the following
	1. `go.db` 
	2. Human and mouse OrgDb packages
	3. `simona` and `simplifyEnrichment` (needed for `ora_simplify`). 

## Running enrichviz

The `enrichviz` command provides a common interface to the ORA and GSEA scripts.

Start a pixi shell, then run commands directly:

```bash
pixi shell
enrichviz --help
```

## Start with differential expression table

If you don't have the functional analysis done, you can start with a differential expression result table.

Enrichviz provides tools to create ORA or GSEA enrichment results.

The command shows to run ORA with gprofiler on a differential expression result table and create a table with  GO functional terms

```bash
enrichviz ora_gprofiler --in data/edger.csv --outdir res --direction yes
```

It produces gprofiler output as a `csv` file in `outdir res`. 

If output file prefix is not specified, the output file will be named as `res/gprofiler_GO.csv`

## Start with ORA/GSEA result table and plot the results.

Enrichviz supports multiple visualization tools. For most of the ORA plotting tools, the only required columns in the input file are `Term` and `Significance`.  `Upset plot` and sematic-space plot (`ssplot`) also need a gene column (eg: geneID,Genes,intersection_genes)

Create a barplot of the top terms

```bash
enrichviz ora_barplot --in res/gprofiler_GO.csv --outdir res --ont BP
```

Get help on an individual command with `-h`.

```bash
enrichviz ora_barplot -h
```

Without a shell, put `--` before the command so pixi does not treat script flags as its own:

```bash
pixi run enrichviz -- ora_barplot -h
```



For plot-specific options and other details, see:

- [ORA visualization guide](docs/ora.md)
- [GSEA visualization guide](docs/gsea.md)



## Available commands

Analysis commands write enrichment tables. Plot commands read those tables and write figures.

### ORA

| Command | Purpose |
|---|---|
| `enrichviz ora_gprofiler` | Run GO enrichment with gProfiler |
| `enrichviz ora_enricher` | Run ORA with clusterProfiler enricher (MSigDB GO or a GMT file) |
| `enrichviz ora_clusterprofiler` | Run GO enrichment with clusterProfiler enrichGO |
| `enrichviz ora_clusterprofiler_kegg` | Run KEGG enrichment with clusterProfiler enrichKEGG |

### ORA plots

| Command | Purpose |
|---|---|
| `enrichviz ora_barplot` | ORA barplot |
| `enrichviz ora_lollipop` | ORA lollipop plot |
| `enrichviz ora_dotplot` | ORA dotplot |
| `enrichviz ora_upset` | UpSet plot |
| `enrichviz ora_ssplot` | Semantic space plot |
| `enrichviz ora_simplify` | Simplify GO terms |

### GSEA

| Command | Purpose |
|---|---|
| `enrichviz gsea_go` | Run GO GSEA |
| `enrichviz gsea_kegg` | Run KEGG GSEA |

### GSEA plots

| Command | Purpose |
|---|---|
| `enrichviz gsea_barplot` | GSEA barplot |
| `enrichviz gsea_lollipop` | GSEA lollipop plot |
| `enrichviz gsea_dotplot` | GSEA dotplot |
| `enrichviz gsea_ridgeplot` | GSEA ridgeplot |
| `enrichviz gsea_volcano` | GSEA volcano plot |
| `enrichviz gsea_nes` | Running ES/NES plot |
| `enrichviz gsea_nes-compare` | Compare NES across GSEA runs |

## ORA or GSEA?

ORA and GSEA answer related but different biological questions.

| | ORA | GSEA |
|---|---|---|
| **Input** | Significant DEGs and a background gene set | All tested genes ranked by a statistic |
| **Question** | Which functions are over-represented among DEGs? | Which functions are enriched toward either end of the ranked gene list? |
| **Direction** | Up- and down-regulated genes can be tested separately | Direction is represented by the ranking statistic and NES |
| **Run** | `enrichviz ora_gprofiler` | `enrichviz gsea_go` or `enrichviz gsea_kegg` |

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


## License

Released under the [MIT License](LICENSE).
