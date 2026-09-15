# Usage:
#   Rscript ora_dotplot.R --in gprofiler_GO.csv --outdir res -n 10 --ont BP
#   Rscript ora_dotplot.R --in kegg.csv --outdir res -n 10 [--prefix NAME]
#
# Dotplot for ORA results from gprofiler, clusterProfiler, fgsea, or a
# similar table (GO, KEGG, or another gene-set database). If a direction
# column is present, top N up and top N down terms are shown together
# (signed gene ratio).

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--in", type = "character", default = "gprofiler_GO.csv", dest = "input_file",
              metavar = "FILE",
              help = "ORA result CSV/TSV or RDS (gprofiler, clusterProfiler, fgsea)"),
  make_option("--outdir", type = "character", default = "res",
              metavar = "DIR", help = "output directory for PDF plots"),
  make_option(c("-n", "--n"), type = "integer", default = 10, dest = "top_n",
              metavar = "N", help = "top terms per direction [default: %default]"),
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
  usage = "\nRscript ora_dotplot.R --in FILE --outdir DIR -n 10 [--ont BP] [--prefix NAME]",
  option_list = option_list,
  description = paste(
    "\nDotplot of top ORA terms from GO, KEGG, or another gene-set database.\n",
    "Accepts gprofiler, clusterProfiler, or fgsea tables (CSV/TSV/RDS).",
    "For GO tables, --ont BP, CC, or MF selects one ontology.",
    "For KEGG and other non-GO tables, ontology filtering is skipped. \n",
    "Up- and down-regulated gene sets are drawn in one plot when direction is present.",
    "--prefix is prepended to output file names.",
    "",
    "Examples:",
    "  Rscript ora_dotplot.R --in gprofiler_GO.csv --outdir res -n 10 --ont BP",
    "  Rscript ora_dotplot.R --in kegg.csv --outdir res -n 10 --prefix KEGG",
    sep = "\n"
  )
)

opt <- parse_args(opt_parser)

if (is.null(opt$input_file) || is.null(opt$outdir)) {
  print_help(opt_parser)
  quit(save = "no", status = 1)
}

input_file <- opt$input_file
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

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

utils_candidates <- c(
  "ora_utils.R",
  file.path("enrich", "ora_utils.R")
)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) == 1) {
  utils_candidates <- c(
    file.path(dirname(normalizePath(sub("^--file=", "", file_arg))), "ora_utils.R"),
    utils_candidates
  )
}
utils_file <- utils_candidates[file.exists(utils_candidates)][1]
if (is.na(utils_file)) {
  stop("Could not find ora_utils.R")
}
source(utils_file)

select_top_terms <- function(df, n) {
  if (!"direction" %in% colnames(df) || all(is.na(df$direction) | df$direction == "")) {
    df$direction <- "none"
  }

  plot_df <- df %>%
    group_by(direction) %>%
    slice_min(order_by = p.adjust, n = n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      score = ifelse(direction == "down", -1, 1) * -log10(pmin(p.adjust, 1)),
      gene_ratio_signed = ifelse(direction == "down", -1, 1) * GeneRatio,
      label = str_to_sentence(Description)
    )

  plot_df <- if (any(is.finite(plot_df$GeneRatio))) {
    arrange(plot_df, GeneRatio, score)
  } else {
    arrange(plot_df, score)
  }

  if (anyDuplicated(plot_df$label) > 0) {
    plot_df <- plot_df %>%
      mutate(label = paste0(label, " (", ID, ")"))
  }

  plot_df %>%
    mutate(label = factor(label, levels = unique(label)))
}

ora_plot_subtitle <- function(plot_df) {
  dirs <- unique(plot_df$direction)
  if (all(c("up", "down") %in% dirs)) {
    "Top positively and negatively enriched terms"
  } else {
    "Top enriched terms"
  }
}

ora_dotplot_plot <- function(plot_df, title) {
  use_ratio <- any(is.finite(plot_df$gene_ratio_signed))
  x_var <- if (use_ratio) "gene_ratio_signed" else "score"
  signed <- ora_has_signed_direction(plot_df)
  x_lab <- if (use_ratio) {
    if (signed) "Signed gene ratio" else "Gene ratio"
  } else if (signed) {
    expression("Signed " * -log[10] * "(p.adjust)")
  } else {
    expression(-log[10] * "(p.adjust)")
  }

  ggplot(plot_df, aes(x = .data[[x_var]], y = label, color = p.adjust, size = Count)) +
    geom_point(alpha = 0.9) +
    ora_zero_vline(plot_df) +
    ora_x_scale(plot_df, plot_df[[x_var]]) +
    scale_color_viridis_c(option = "magma", direction = -1, trans = "log10") +
    labs(
      title = title,
      subtitle = ora_plot_subtitle(plot_df),
      x = x_lab,
      y = NULL,
      color = "p.adjust",
      size = "Gene count"
    ) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 10))
}

message("Input:  ", input_file)
message("Output: ", RESULTS_DIR)
message("Prefix: ", if (nzchar(prefix)) prefix else "(none)")
message("N:      ", top_n)

loaded <- prepare_ora_for_plot(input_file, ont = ONT, fdr_cutoff = fdr_cutoff)
plot_meta <- loaded$meta
ont_filter <- loaded$ont_filter
ora_df <- loaded$df

message("Source: ", plot_meta$source)
message("ONT:    ", if (length(ont_filter) == 1) ont_filter else plot_meta$source)

if (nrow(ora_df) == 0) {
  stop(
    "No ", plot_meta$term_label, " terms with FDR < ", fdr_cutoff, " in ", input_file,
    call. = FALSE
  )
}

plot_df <- select_top_terms(ora_df, top_n)
if (nrow(plot_df) == 0) {
  stop("No terms left to plot after selecting the top ", top_n, ".", call. = FALSE)
}

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

plot_csv <- file.path(RESULTS_DIR, ora_out_name(paste0(plot_meta$tag, "_dotplot_table.csv"), prefix))
write_csv(plot_df %>% mutate(label = as.character(label)), plot_csv)

n_rows <- nrow(plot_df)
plot_height <- max(5, 0.27 * n_rows + 2)

dot_file <- file.path(RESULTS_DIR, ora_out_name(paste0(plot_meta$tag, "_dotplot.pdf"), prefix))

ggsave(
  dot_file,
  plot = ora_dotplot_plot(plot_df, plot_meta$title),
  width = 9,
  height = plot_height
)

dir_counts <- table(plot_df$direction)

message("================")
message("Significant ", plot_meta$term_label, " terms (FDR < ", fdr_cutoff, "): ", nrow(ora_df))
message(
  "Plotted: ",
  paste(paste0(names(dir_counts), "=", as.integer(dir_counts)), collapse = ", ")
)
message("CSV files written:")
message("  ", plot_csv)
message("Plots created:")
message("Dotplot of significant ", plot_meta$term_label, " terms: ", dot_file)
