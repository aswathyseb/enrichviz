# enrichviz

R scripts for running and plotting gene-set enrichment from differential-expression results. Two complementary analyses are supported:

- **Over-representation analysis (ORA)** tests whether annotated functions occur more often among significant genes than expected by chance.
- **Gene set enrichment analysis (GSEA)** tests whether genes in a function concentrate toward the top or bottom of a ranked list of all tested genes.

Both tracks accept GO and KEGG results and write PDF plots plus the tables used to draw them. Plot-level details, options, and examples live in the analysis-specific guides:

- [Visualization of enriched terms from ORA](ora/README_ORA.md)
- [Visualization of enriched terms from GSEA](gsea/README_GSEA.md)

Install the R environment with [install.md](install.md). After that, run scripts from the repository root, either inside `pixi shell` or with `pixi run Rscript …`.

```bash
pixi run Rscript ora/ora_barplot.R -h
pixi run Rscript gsea/gsea_barplot.R -h
```

## ORA or GSEA

| | ORA | GSEA |
|---|---|---|
| Input genes | Significant DEGs (and a background) | All tested genes, ranked |
| Typical question | Which functions are over-represented among DEGs? | Do function members pile up among the strongest up- or down-regulated genes? |
| Direction | Optional: test up and down DEGs separately | Encoded in the ranking statistic and NES |
| Run the analysis | `gprof_GO.R` or `clusterProf_GO.R` | `gsea_GO.R` / `gsea_KEGG.R` (fgsea) or the clusterProfiler wrappers |

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
