# Usage:
#   Rscript gsea_ridgeplot.R --in gsea_out/fgsea_GO_all.csv --rank gsea_out/ranked_genes_GO.csv --outdir gsea_out -n 10 --ont BP [--prefix NAME]
#   Rscript gsea_ridgeplot.R --in t5/fgsea_KEGG.csv --rank t5/ranked_genes_KEGG.csv --outdir gsea_out -n 10 [--prefix kegg]
#
# Ridgeplot for GO or KEGG GSEA results from fgsea or clusterProfiler.
# Significant terms (FDR < 0.05) are split by NES sign; top N up and
# top N down terms are shown together. Each ridge is the rank-score
# density of that term's leading-edge genes.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--in", type = "character", default = "gsea_out/fgsea_GO_all.csv",
              dest = "input_file", metavar = "FILE",
              help = "GSEA result CSV/TSV or RDS (fgsea or clusterProfiler)"),
  make_option("--rank", type = "character", default = "gsea_out/ranked_genes_GO.csv",
              dest = "rank_file", metavar = "FILE",
              help = "ranked gene CSV/TSV or RDS (rank_score + gene_symbol/ensembl)"),
  make_option("--outdir", type = "character", default = "gsea_out",
              metavar = "DIR", help = "output directory for PDF plots"),
  make_option(c("-n", "--n"), type = "integer", default = 10, dest = "top_n",
              metavar = "N", help = "top terms per NES direction [default: %default]"),
  make_option("--ont", type = "character", default = "BP",
              metavar = "BP|CC|MF|KEGG|all",
              help = paste(
                "GO ontology to plot, or KEGG/all to skip GO filtering",
                "[default: %default]"
              )),
  make_option("--prefix", type = "character", default = "",
              metavar = "NAME",
              help = "optional prefix for output file names [default: none]")
)

opt_parser <- OptionParser(
  usage = "\nRscript gsea_ridgeplot.R --in FILE --rank FILE --outdir DIR -n 10 [--ont BP] [--prefix NAME]\n",
  option_list = option_list,
  description = paste(
    "Ridgeplot of top GSEA terms from GO or KEGG.",
    "Accepts fgsea or clusterProfiler tables (CSV/TSV/RDS).\n",
    "For GO tables, --ont BP, CC, or MF selects one ontology.",
    "For KEGG tables, ontology filtering is skipped.\n",
    "Each ridge shows the ranking-statistic density of leading-edge genes.",
    "Positive- and negative-NES gene sets are drawn on one plot.",
    "",
    "Examples:",
    "  Rscript gsea_ridgeplot.R --in gsea_out/fgsea_GO_all.csv --rank gsea_out/ranked_genes_GO.csv --outdir gsea_out -n 10 --ont BP",
    "  Rscript gsea_ridgeplot.R --in t5/fgsea_KEGG.csv --rank t5/ranked_genes_KEGG.csv --outdir gsea_out -n 10 --prefix kegg",
    sep = "\n"
  )
)

opt <- parse_args(opt_parser)

if (is.null(opt$input_file) || is.null(opt$outdir)) {
  print_help(opt_parser)
  quit(save = "no", status = 1)
}

input_file <- opt$input_file
rank_file <- opt$rank_file
RESULTS_DIR <- opt$outdir
top_n <- as.integer(opt$top_n)
ONT <- toupper(opt$ont)
prefix <- if (is.null(opt$prefix)) "" else as.character(opt$prefix)
fdr_cutoff <- 0.05

if (!ONT %in% c("BP", "CC", "MF", "KEGG", "ALL")) {
  stop("--ont must be BP, CC, MF, KEGG, or all", call. = FALSE)
}
if (!is.finite(top_n) || top_n < 1) {
  stop("-n must be a positive integer", call. = FALSE)
}
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}
if (is.null(rank_file) || !file.exists(rank_file)) {
  stop("Ranked-gene file not found: ", rank_file, call. = FALSE)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(ggridges)
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

select_top_terms <- function(df, n) {
  plot_df <- df %>%
    group_by(direction) %>%
    slice_max(order_by = abs(NES), n = n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(label = str_to_sentence(Description)) %>%
    arrange(NES)

  if (anyDuplicated(plot_df$label) > 0) {
    plot_df <- plot_df %>%
      mutate(label = paste0(label, " (", ID, ")"))
  }

  plot_df %>%
    mutate(label = factor(label, levels = unique(label)))
}

gsea_plot_subtitle <- function(plot_df) {
  dirs <- unique(plot_df$direction)
  if (all(c("up", "down") %in% dirs)) {
    "Leading-edge rank scores for top positively and negatively enriched terms"
  } else if ("up" %in% dirs) {
    "Leading-edge rank scores for top positively enriched terms"
  } else if ("down" %in% dirs) {
    "Leading-edge rank scores for top negatively enriched terms"
  } else {
    "Leading-edge rank scores for top enriched terms"
  }
}

gsea_ridgeplot_plot <- function(ridge_df, title, subtitle) {
  ggplot(ridge_df, aes(x = rank_score, y = label, fill = p.adjust)) +
    geom_density_ridges(
      scale = 1.1,
      alpha = 0.85,
      rel_min_height = 0.01,
      color = "grey30",
      linewidth = 0.2
    ) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    scale_fill_viridis_c(option = "magma", direction = -1, trans = "log10") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Rank score: sign(log2FC) * -log10(PValue)",
      y = NULL,
      fill = "p.adjust"
    ) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 10))
}

message("Input:  ", input_file)
message("Ranks:  ", rank_file)
message("Output: ", RESULTS_DIR)
message("Prefix: ", if (nzchar(prefix)) prefix else "(none)")
message("N:      ", top_n)

loaded <- prepare_gsea_for_plot(input_file, ont = ONT, fdr_cutoff = fdr_cutoff)
plot_meta <- loaded$meta
ont_filter <- loaded$ont_filter
gsea_df <- loaded$df

message("Source: ", plot_meta$source)
message("ONT:    ", if (length(ont_filter) == 1) ont_filter else plot_meta$source)
if (nrow(gsea_df) == 0) {
  stop(
    "No ", plot_meta$term_label, " terms with FDR < ", fdr_cutoff, " in ", input_file,
    call. = FALSE
  )
}

plot_df <- select_top_terms(gsea_df, top_n)
if (nrow(plot_df) == 0) {
  stop("No terms left to plot after selecting the top ", top_n, ".", call. = FALSE)
}

ranks <- read_gsea_rank_table(rank_file)
ridge_df <- expand_gsea_leading_edge(plot_df, ranks)
if (nrow(ridge_df) == 0) {
  stop(
    "No leading-edge genes overlapped the ranked list in ", rank_file,
    ". Use ranked_genes_GO.csv or ranked_genes_KEGG.csv so symbols and Ensembl IDs can both be matched.",
    call. = FALSE
  )
}

n_genes <- ridge_df %>% count(label, name = "n_matched")
keep_labels <- n_genes$label[n_genes$n_matched >= 3]
if (length(keep_labels) == 0) {
  stop("No selected terms had at least 3 leading-edge genes with rank scores.", call. = FALSE)
}
dropped <- setdiff(levels(plot_df$label), as.character(keep_labels))
if (length(dropped) > 0) {
  message(
    "Dropped ", length(dropped),
    " term(s) with fewer than 3 matched leading-edge genes."
  )
}

ridge_df <- ridge_df %>%
  filter(label %in% keep_labels) %>%
  mutate(label = factor(as.character(label), levels = levels(plot_df$label)[levels(plot_df$label) %in% keep_labels]))

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

plot_csv <- file.path(RESULTS_DIR, gsea_out_name(paste0(plot_meta$tag, "_ridgeplot_table.csv"), prefix))
write_csv(ridge_df %>% mutate(label = as.character(label)), plot_csv)

n_rows <- nlevels(ridge_df$label)
plot_height <- max(6, 0.42 * n_rows + 2)

ridge_file <- file.path(RESULTS_DIR, gsea_out_name(paste0(plot_meta$tag, "_ridgeplot.pdf"), prefix))

ggsave(
  ridge_file,
  plot = gsea_ridgeplot_plot(ridge_df, plot_meta$title, gsea_plot_subtitle(plot_df)),
  width = 9,
  height = plot_height
)

dir_counts <- table(plot_df$direction[plot_df$label %in% levels(ridge_df$label)])

message("================")
message("Significant ", plot_meta$term_label, " terms (FDR < ", fdr_cutoff, "): ", nrow(gsea_df))
message(
  "Plotted: ",
  paste(paste0(names(dir_counts), "=", as.integer(dir_counts)), collapse = ", ")
)
message("CSV files written:")
message("  ", plot_csv)
message("Plots created:")
message("Ridgeplot of significant ", plot_meta$term_label, " terms: ", ridge_file)
