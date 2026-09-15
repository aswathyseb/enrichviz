# Usage:
#   Rscript gsea_dotplot.R --in gsea_out/fgsea_GO_all.csv --outdir gsea_out -n 10 --ont BP [--prefix NAME]
#   Rscript gsea_dotplot.R --in t5/fgsea_KEGG.csv --outdir gsea_out -n 10 [--prefix kegg]
#
# Dotplot for GO or KEGG GSEA results from fgsea or clusterProfiler.
# Significant terms (FDR < 0.05) are split by NES sign; top N up and
# top N down terms are shown together (x = signed leading-edge ratio).

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--in", type = "character", default = "gsea_out/fgsea_GO_all.csv",
              dest = "input_file", metavar = "FILE",
              help = "GSEA result CSV/TSV or RDS (fgsea or clusterProfiler)"),
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
  usage = "\nRscript gsea_dotplot.R --in FILE --outdir DIR -n 10 [--ont BP] [--prefix NAME]\n",
  option_list = option_list,
  description = paste(
    "Dotplot of top GSEA terms from GO or KEGG.",
    "Accepts fgsea or clusterProfiler tables (CSV/TSV/RDS).\n",
    "For GO tables, --ont BP, CC, or MF selects one ontology.",
    "For KEGG tables, ontology filtering is skipped.\n",
    "x is signed leading-edge ratio (leading-edge genes / set size).",
    "Positive- and negative-NES gene sets are drawn on one plot.",
    "",
    "Examples:",
    "  Rscript gsea_dotplot.R --in gsea_out/fgsea_GO_all.csv --outdir gsea_out -n 10 --ont BP",
    "  Rscript gsea_dotplot.R --in t5/fgsea_KEGG.csv --outdir gsea_out -n 10 --prefix kegg",
    "  Rscript gsea_dotplot.R --in t5/gseaKEGG_clusterprof.csv --outdir gsea_out -n 10 --prefix kegg",
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
    mutate(
      leadingEdgeRatio_signed = ifelse(direction == "down", -1, 1) * leadingEdgeRatio,
      label = str_to_sentence(Description)
    )

  plot_df <- if (any(is.finite(plot_df$leadingEdgeRatio_signed))) {
    arrange(plot_df, leadingEdgeRatio_signed, NES)
  } else {
    arrange(plot_df, NES)
  }

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
    "Top positively and negatively enriched terms"
  } else if ("up" %in% dirs) {
    "Top positively enriched terms"
  } else if ("down" %in% dirs) {
    "Top negatively enriched terms"
  } else {
    "Top enriched terms"
  }
}

gsea_dotplot_plot <- function(plot_df, title) {
  use_ratio <- any(is.finite(plot_df$leadingEdgeRatio_signed))
  x_var <- if (use_ratio) "leadingEdgeRatio_signed" else "NES"
  x_lab <- if (use_ratio) {
    "Signed leading-edge ratio"
  } else {
    "Normalized enrichment score (NES)"
  }

  ggplot(plot_df, aes(x = .data[[x_var]], y = label, color = p.adjust, size = leadingEdgeCount)) +
    geom_point(alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    scale_color_viridis_c(option = "magma", direction = -1, trans = "log10") +
    labs(
      title = title,
      subtitle = gsea_plot_subtitle(plot_df),
      x = x_lab,
      y = NULL,
      color = "p.adjust",
      size = "Leading-edge genes"
    ) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 10))
}

message("Input:  ", input_file)
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
if (!any(is.finite(plot_df$leadingEdgeRatio))) {
  message("Leading-edge ratio is missing; falling back to NES on the x-axis.")
}

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

plot_csv <- file.path(RESULTS_DIR, gsea_out_name(paste0(plot_meta$tag, "_dotplot_table.csv"), prefix))
write_csv(plot_df %>% mutate(label = as.character(label)), plot_csv)

n_rows <- nrow(plot_df)
plot_height <- max(5, 0.27 * n_rows + 2)

dot_file <- file.path(RESULTS_DIR, gsea_out_name(paste0(plot_meta$tag, "_dotplot.pdf"), prefix))

ggsave(
  dot_file,
  plot = gsea_dotplot_plot(plot_df, plot_meta$title),
  width = 9,
  height = plot_height
)

dir_counts <- table(plot_df$direction)

message("================")
message("Significant ", plot_meta$term_label, " terms (FDR < ", fdr_cutoff, "): ", nrow(gsea_df))
message(
  "Plotted: ",
  paste(paste0(names(dir_counts), "=", as.integer(dir_counts)), collapse = ", ")
)
message("CSV files written:")
message("  ", plot_csv)
message("Plots created:")
message("Dotplot of significant ", plot_meta$term_label, " terms: ", dot_file)
