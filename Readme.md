# enrichviz

R scripts for running and plotting gene-set enrichment from differential-expression results. Two complementary analyses are supported:

- **Over-representation analysis (ORA)** tests whether annotated functions occur more often among significant genes than expected by chance.
- **Gene set enrichment analysis (GSEA)** tests whether genes in a function concentrate toward the top or bottom of a ranked list of all tested genes.

Both tracks accept GO and KEGG results and write PDF plots plus the tables used to draw them. Plot-level details, options, and examples live in the analysis-specific guides:

- [Visualization of enriched terms from ORA](ora/README_ORA.md)
- [Visualization of enriched terms from GSEA](gsea/README_GSEA.md)

## Install

The environment is defined by `pixi.toml` and pinned in `pixi.lock`. You need [pixi](https://pixi.sh/).

```bash
git clone https://github.com/aswathyseb/enrichviz.git
cd enrichviz
pixi install
pixi run setup-annot
```

`setup-annot` downloads `go.db` and the human/mouse OrgDb packages. Pixi does not run conda post-link scripts, so `pixi install` alone leaves those stubs empty.

`ora simplify` also needs `simona` and `simplifyEnrichment`. Those steps, and rebuilding the environment from scratch, are in [install.md](install.md).

## Run

`enrichviz` maps subcommands to the R scripts. Script flags (`--in`, `--outdir`, `-h`, …) are unchanged.

Inside `pixi shell`:

```bash
pixi shell
enrichviz --help
enrichviz ora barplot -h
enrichviz ora barplot --in ora_out/gprofiler_GO.csv --outdir ora_out -n 10 --ont BP
enrichviz gsea go --in edger.csv --outdir gsea_out --ont all
```

Inside `pixi shell`, Tab completes tracks and subcommands (`enrichviz ora gpr<Tab>` → `gprofiler`). After the subcommand, Tab completes file paths. Start a new `pixi shell` after pulling this change.

Without a shell, start the same command with `pixi run enrichviz --` so pixi does not treat `--in` as its own flag:

```bash
pixi run enrichviz -- ora barplot -h
pixi run enrichviz -- gsea barplot --in gsea_out/fgsea_GO_all.csv --outdir gsea_out -n 10 --ont BP
```

| Command | Script |
|---|---|
| `enrichviz ora gprofiler` | `ora/gprof_GO.R` |
| `enrichviz ora barplot` | `ora/ora_barplot.R` |
| `enrichviz ora lollipop` | `ora/ora_lollipop.R` |
| `enrichviz ora dotplot` | `ora/ora_dotplot.R` |
| `enrichviz ora upset` | `ora/ora_upsetplot.R` |
| `enrichviz ora ssplot` | `ora/ora_ssplot.R` |
| `enrichviz ora simplify` | `ora/simplify_GO.R` |
| `enrichviz gsea go` | `gsea/gsea_GO.R` |
| `enrichviz gsea kegg` | `gsea/gsea_KEGG.R` |
| `enrichviz gsea barplot` | `gsea/gsea_barplot.R` |
| `enrichviz gsea lollipop` | `gsea/gsea_lollipop.R` |
| `enrichviz gsea dotplot` | `gsea/gsea_dotplot.R` |
| `enrichviz gsea ridgeplot` | `gsea/gsea_ridgeplot.R` |
| `enrichviz gsea volcano` | `gsea/gsea_volcano.R` |
| `enrichviz gsea nes` | `gsea/gsea_NES.R` |
| `enrichviz gsea nes-compare` | `gsea/gsea_NES_compare.R` |

## ORA or GSEA

| | ORA | GSEA |
|---|---|---|
| Input genes | Significant DEGs (and a background) | All tested genes, ranked |
| Typical question | Which functions are over-represented among DEGs? | Do function members pile up among the strongest up- or down-regulated genes? |
| Direction | Optional: test up and down DEGs separately | Encoded in the ranking statistic and NES |
| Run the analysis | `enrichviz ora gprofiler` | `enrichviz gsea go` / `enrichviz gsea kegg` |

Use ORA when you already have a DEG list and want functions that are over-represented in that list. Use GSEA when you want to use the full ranked list, including genes that do not pass a significance cutoff.

Plotting scripts read CSV, TSV, or RDS tables from gProfiler, clusterProfiler, or fgsea.

## Choosing a visualization

Several plot types appear in both tracks. They answer the same kind of question, but the axes differ: ORA emphasizes significance and gene overlap among DEGs, while GSEA emphasizes NES, leading-edge genes, and position in the ranked list.

**Overview of the strongest terms**

- **Barplot** — magnitude and direction of the top terms.
- **Lollipop** — the same ranking, with point size for gene-set or gene count.
- **Dotplot** — fraction of the query (ORA gene ratio) or of the gene set that drives enrichment (GSEA leading-edge ratio).

**How terms relate to each other** (ORA)

- **UpSet plot** — shared query genes among selected terms.
- **Semantic space plot** — terms arranged by Jaccard overlap of query genes.
- **simplifyGO** — GO terms grouped by semantic similarity in the ontology.

**Where the signal sits in the ranked list** (GSEA)

- **Ridgeplot** — distribution of leading-edge ranking scores.
- **Volcano** — NES versus adjusted *P* for all tested terms.
- **Running ES / NES** — enrichment trajectory for one term.
- **NES comparison** — the same term across multiple GSEA runs.

For term-selection rules, plot elements, and command lines, see [ora/README_ORA.md](ora/README_ORA.md) and [gsea/README_GSEA.md](gsea/README_GSEA.md).

Released under the [MIT License](LICENSE).

