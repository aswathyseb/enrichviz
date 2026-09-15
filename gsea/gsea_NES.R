# Usage:
#   Rscript gsea_NES.R --in gsea_out/fgsea_GO_all.rds --GO GO:0051607 --outdir gsea_out [--prefix NAME]
#   Rscript gsea_NES.R --in t5/fgsea_KEGG.rds --term hsa03010 --outdir gsea_out [--prefix kegg]
#
# Running-enrichment plot for one GO or KEGG term (gseaNb-style).
# With --expr (DE table containing normalized counts after falsePos),
# a z-scored leading-edge heatmap is stacked under the NES curve.
# Without --expr, only the running ES and ranked-list panels are drawn.



suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--in", type = "character", default = "gsea_out/fgsea_GO_all.rds",
              dest = "input_file", metavar = "FILE",
              help = "GSEA result RDS or CSV (fgsea or clusterProfiler)"),
  make_option("--GO", type = "character", default = NULL, dest = "GO_ID",
              metavar = "ID", help = "GO or KEGG term to plot (required unless --term)"),
  make_option("--term", type = "character", default = NULL, dest = "term_id",
              metavar = "ID", help = "alias for --GO"),
  make_option("--expr", type = "character", default = NULL, dest = "expr_file",
              metavar = "FILE",
              help = "optional DE CSV with normalized counts after falsePos"),
  make_option("--rank", type = "character", default = "gsea_out/ranked_genes_GO.csv",
              dest = "rank_file", metavar = "FILE",
              help = "ranked-gene file used if the GSEA RDS has no stats"),
  make_option("--outdir", type = "character", default = "gsea_out",
              metavar = "DIR", help = "output directory for PDF plots"),
  make_option("--prefix", type = "character", default = "",
              metavar = "NAME",
              help = "optional prefix for output file names [default: none]")
)

opt_parser <- OptionParser(
  usage = "\nRscript gsea_NES.R --in FILE --GO ID [--expr FILE] --outdir DIR [--prefix NAME]\n",
  option_list = option_list,
  description = paste(
    "Running enrichment (NES) plot for one GO or KEGG term.",
    "Prefers an fgsea/clusterProfiler RDS so the full gene set can be used.",
    "If --expr is given, a z-scored heatmap of leading-edge genes is added.",
    "",
    "Examples:",
    "  Rscript gsea_NES.R --in gsea_out/fgsea_GO_all.csv --GO GO:0051607 --outdir gsea_out",
    "  Rscript gsea_NES.R --in t5/fgsea_KEGG.rds --term hsa03010 --outdir gsea_out --prefix kegg",
    sep = "\n"
  )
)

opt <- parse_args(opt_parser)

term_id <- opt$term_id
if (is.null(term_id) || !nzchar(term_id)) {
  term_id <- opt$GO_ID
}

if (is.null(opt$input_file) || is.null(opt$outdir) || is.null(term_id) || !nzchar(term_id)) {
  print_help(opt_parser)
  quit(save = "no", status = 1)
}

input_file <- opt$input_file
GO_ID <- term_id
expr_file <- opt$expr_file
rank_file <- opt$rank_file
RESULTS_DIR <- opt$outdir
prefix <- if (is.null(opt$prefix)) "" else as.character(opt$prefix)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}
if (!is.null(expr_file) && !nzchar(expr_file)) {
  expr_file <- NULL
}
if (!is.null(expr_file) && !file.exists(expr_file)) {
  stop("Expression file not found: ", expr_file, call. = FALSE)
}
if (!is.null(rank_file) && !file.exists(rank_file)) {
  rank_file <- NULL
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(patchwork)
})

utils_candidates <- c(
  "gsea_utils.R",
  file.path("enrich", "gsea_utils.R")
)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) == 1) {
  utils_candidates <- c(
    file.path(dirname(normalizePath(sub("^--file=", "", file_arg))), "gsea_utils.R"),
    utils_candidates
  )
}
utils_file <- utils_candidates[file.exists(utils_candidates)][1]
if (is.na(utils_file)) {
  stop("Could not find gsea_utils.R")
}
source(utils_file)

gsea_es_plot <- function(es_df, term) {
  title <- str_to_title(term$Description[[1]])
  nes <- term$NES[[1]]
  fdr <- term$p.adjust[[1]]
  stats_lab <- sprintf("NES = %.2f   FDR = %s", nes, format(fdr, digits = 2, scientific = TRUE))

  ggplot(es_df, aes(x = rank, y = running_es)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, color = "grey40") +
    geom_line(color = "#4DAF4A", linewidth = 0.7) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, max(es_df$rank))) +
    labs(
      title = title,
      subtitle = stats_lab,
      x = NULL,
      y = "Running Enrichment Score"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

gsea_rank_plot <- function(es_df) {
  ticks <- es_df %>% filter(in_set)
  ggplot(es_df, aes(x = rank, y = 1, fill = stat)) +
    geom_raster() +
    geom_segment(
      data = ticks,
      aes(x = rank, xend = rank, y = 0.55, yend = 1.45),
      inherit.aes = FALSE,
      color = "grey15",
      linewidth = 0.15
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max(es_df$rank)), expand = FALSE) +
    labs(x = "Ranked List", y = NULL, fill = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none"
    )
}

heatmap_from_expr <- function(expr_file, genes, stats) {
  de <- read_csv(expr_file, show_col_types = FALSE)
  expr <- parse_expr(de)
  sample_cols <- colnames(expr)[-1]
  if (length(sample_cols) == 0) {
    stop("No sample columns were found after 'falsePos' in ", expr_file, call. = FALSE)
  }

  expr_key <- expr %>%
    mutate(gene_key = toupper(gene_symbol)) %>%
    distinct(gene_key, .keep_all = TRUE)

  ens_map <- NULL
  if (all(c("name", "gene") %in% colnames(de))) {
    ens_map <- de %>%
      transmute(
        ensembl_gene = sub("\\..*$", "", as.character(name)),
        gene_symbol = as.character(gene)
      ) %>%
      filter(!is.na(gene_symbol), gene_symbol != "") %>%
      distinct(ensembl_gene, .keep_all = TRUE)
  }

  genes <- genes[nzchar(genes) & !is.na(genes)]
  genes <- genes[!duplicated(toupper(genes))]
  matched <- tibble::tibble(query = genes, gene_key = toupper(genes)) %>%
    left_join(expr_key %>% select(gene_key, gene_symbol), by = "gene_key")

  if (!is.null(ens_map) && any(is.na(matched$gene_symbol))) {
    by_ens <- matched %>%
      filter(is.na(gene_symbol)) %>%
      select(query, gene_key) %>%
      left_join(
        ens_map %>% transmute(gene_key = toupper(ensembl_gene), gene_symbol),
        by = "gene_key"
      )
    matched$gene_symbol[is.na(matched$gene_symbol)] <- by_ens$gene_symbol[
      match(matched$query[is.na(matched$gene_symbol)], by_ens$query)
    ]
  }

  matched <- matched %>%
    filter(!is.na(gene_symbol), gene_symbol %in% expr$gene_symbol) %>%
    distinct(gene_symbol, .keep_all = TRUE)

  if (nrow(matched) == 0) {
    return(NULL)
  }

  matched$rank_in_list <- vapply(seq_len(nrow(matched)), function(i) {
    keys <- toupper(c(matched$query[[i]], matched$gene_symbol[[i]]))
    hits <- match(keys, toupper(names(stats)))
    hits <- hits[!is.na(hits)]
    if (length(hits) == 0) Inf else min(hits)
  }, numeric(1))
  matched <- matched %>% arrange(rank_in_list, gene_symbol)

  mat <- expr[match(matched$gene_symbol, expr$gene_symbol), sample_cols, drop = FALSE]
  mat <- as.matrix(mat)
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0

  heat_df <- as.data.frame(z, stringsAsFactors = FALSE)
  heat_df$gene <- factor(matched$gene_symbol, levels = matched$gene_symbol)
  heat_df %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "sample",
      values_to = "zscore"
    ) %>%
    mutate(sample = factor(sample, levels = rev(sample_cols)))
}

gsea_heatmap_plot <- function(heat_df) {
  ggplot(heat_df, aes(x = gene, y = sample, fill = zscore)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Z-Score"
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      axis.text.y = element_text(size = 8),
      axis.ticks = element_blank(),
      legend.key.height = unit(0.45, "cm"),
      plot.margin = margin(4, 8, 4, 4)
    )
}

message("Input:  ", input_file)
message("Term:   ", GO_ID)
message("Expr:   ", if (is.null(expr_file)) "(none)" else expr_file)
message("Output: ", RESULTS_DIR)
if (nzchar(prefix)) message("Prefix: ", prefix)

nes_in <- load_gsea_nes_term(input_file, GO_ID, rank_file = rank_file)
if (is.null(nes_in$stats)) {
  stop(
    "Could not find a ranked gene list. Pass an fgsea/clusterProfiler RDS ",
    "or --rank ranked_genes_GO.csv / ranked_genes_KEGG.csv.",
    call. = FALSE
  )
}

es_df <- gsea_running_es(nes_in$stats, nes_in$pathway_genes)
term <- nes_in$term

heat_df <- NULL
if (!is.null(expr_file)) {
  ht_genes <- if (length(nes_in$leading_edge) > 0) nes_in$leading_edge else nes_in$pathway_genes
  heat_df <- heatmap_from_expr(expr_file, ht_genes, nes_in$stats)
  if (is.null(heat_df) || nrow(heat_df) == 0) {
    message("No leading-edge genes overlapped the expression table; drawing the NES plot only.")
    heat_df <- NULL
  }
}

PLOT_WIDTH <- 8.5
PLOT_HEIGHT <- 5.5

es_p <- gsea_es_plot(es_df, term)
rank_p <- gsea_rank_plot(es_df)

if (is.null(heat_df)) {
  combined <- es_p / rank_p + plot_layout(heights = c(3.2, 0.7))
  plot_height <- PLOT_HEIGHT
  plot_width <- PLOT_WIDTH
} else {
  n_genes <- nlevels(heat_df$gene)
  n_samples <- nlevels(heat_df$sample)
  ht_height <- max(2.2, min(3.6, 0.18 * n_samples + 1.0))
  combined <- es_p / rank_p / gsea_heatmap_plot(heat_df) +
    plot_layout(heights = c(2.5, 0.55, ht_height))
  plot_height <- 4.0 + ht_height
  plot_width <- max(8.5, min(14, 0.16 * n_genes + 3.5))
}

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

go_tag <- gsea_term_tag(nes_in$go_id)
plot_csv <- file.path(RESULTS_DIR, gsea_out_name(paste0("gsea_", go_tag, "_NES_table.csv"), prefix))
write_csv(es_df %>% filter(in_set), plot_csv)

nes_file <- file.path(RESULTS_DIR, gsea_out_name(paste0("gsea_", go_tag, "_NES.pdf"), prefix))
ggsave(nes_file, plot = combined, width = plot_width, height = plot_height)

message("================")
message("Term: ", term$ID[[1]], " — ", term$Description[[1]])
message("NES:  ", signif(term$NES[[1]], 4), "   FDR: ", signif(term$p.adjust[[1]], 3))
message("Gene-set genes in ranked list: ", sum(es_df$in_set), " / ", nrow(es_df))
if (!is.null(heat_df)) {
  message("Heatmap genes: ", nlevels(heat_df$gene))
}
message("CSV files written:")
message("  ", plot_csv)
message("Plots created:")
message("NES plot: ", nes_file)
