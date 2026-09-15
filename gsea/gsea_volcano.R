# Usage:
#   Rscript gsea_volcano.R --in gsea_out/fgsea_GO_all.csv --outdir gsea_out -n 10 --ont BP [--prefix NAME]
#   Rscript gsea_volcano.R --in t5/fgsea_KEGG.csv --outdir gsea_out -n 10 [--prefix kegg]
#
# Volcano plot for GO or KEGG GSEA results from fgsea or clusterProfiler.
# All tested terms in the ontology (or all KEGG pathways) are drawn
# (x = NES, y = -log10 FDR). Top N significant terms per NES direction
# are labeled.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--in", type = "character", default = "gsea_out/fgsea_GO_all.csv",
              dest = "input_file", metavar = "FILE",
              help = "GSEA result CSV/TSV or RDS (fgsea or clusterProfiler)"),
  make_option("--outdir", type = "character", default = "gsea_out",
              metavar = "DIR", help = "output directory for PDF plots"),
  make_option(c("-n", "--n"), type = "integer", default = 10, dest = "top_n",
              metavar = "N", help = "terms to label per NES direction [default: %default]"),
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
  usage = "\nRscript gsea_volcano.R --in FILE --outdir DIR -n 10 [--ont BP] [--prefix NAME]\n",
  option_list = option_list,
  description = paste(
    "Volcano plot of GSEA terms from GO or KEGG.",
    "Accepts fgsea or clusterProfiler tables (CSV/TSV/RDS).\n",
    "For GO tables, --ont BP, CC, or MF selects one ontology.",
    "For KEGG tables, ontology filtering is skipped.",
    "Plots all tested terms; labels the top N significant terms per NES direction.",
    "",
    "Examples:",
    "  Rscript gsea_volcano.R --in gsea_out/fgsea_GO_all.csv --outdir gsea_out -n 10 --ont BP",
    "  Rscript gsea_volcano.R --in t5/fgsea_KEGG.csv --outdir gsea_out -n 10 --prefix kegg",
    "  Rscript gsea_volcano.R --in t5/gseaKEGG_clusterprof.csv --outdir gsea_out -n 10 --prefix kegg",
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
  library(ggrepel)
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

prepare_volcano_table <- function(path, ont, fdr_cutoff) {
  loaded <- prepare_gsea_for_plot(path, ont = ont, fdr_cutoff = Inf)
  df <- loaded$df %>%
    mutate(
      direction = ifelse(NES > 0, "up", "down"),
      significant = p.adjust < fdr_cutoff,
      status = case_when(
        significant & direction == "up" ~ "up",
        significant & direction == "down" ~ "down",
        TRUE ~ "ns"
      ),
      neglog10_fdr = -log10(pmin(p.adjust, 1)),
      label = str_to_sentence(Description)
    )
  list(df = df, meta = loaded$meta, ont_filter = loaded$ont_filter)
}

select_labels <- function(df, n) {
  labeled <- df %>%
    filter(significant) %>%
    group_by(direction) %>%
    slice_max(order_by = abs(NES), n = n, with_ties = FALSE) %>%
    ungroup()

  if (anyDuplicated(labeled$label) > 0) {
    labeled <- labeled %>%
      mutate(label = paste0(label, " (", ID, ")"))
  }
  labeled
}

gsea_volcano_plot <- function(plot_df, label_df, title, fdr_cutoff) {
  status_colors <- c(
    down = "#327eba",
    ns = "grey75",
    up = "#e06663"
  )

  ggplot(plot_df, aes(x = NES, y = neglog10_fdr)) +
    geom_point(aes(color = status), alpha = 0.75, size = 1.6) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2, linewidth = 0.4) +
    geom_text_repel(
      data = label_df,
      aes(label = label, color = status),
      size = 3,
      max.overlaps = 40,
      min.segment.length = 0,
      segment.color = "grey40",
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = status_colors,
      breaks = c("up", "down", "ns"),
      labels = c("Up (p.adjust < 0.05)", "Down (p.adjust < 0.05)", "Not significant")
    ) +
    labs(
      title = title,
      subtitle = "All tested terms; labels are the top NES terms per direction",
      x = "Normalized enrichment score (NES)",
      y = expression(-log[10] * "(p.adjust)"),
      color = NULL
    ) +
    theme_classic(base_size = 11)
}

message("Input:  ", input_file)
message("Output: ", RESULTS_DIR)
message("Prefix: ", if (nzchar(prefix)) prefix else "(none)")
message("N:      ", top_n)

loaded <- prepare_volcano_table(input_file, ont = ONT, fdr_cutoff = fdr_cutoff)
plot_df <- loaded$df
plot_meta <- loaded$meta
ont_filter <- loaded$ont_filter

message("Source: ", plot_meta$source)
message("ONT:    ", if (length(ont_filter) == 1) ont_filter else plot_meta$source)
if (nrow(plot_df) == 0) {
  stop(
    "No ", plot_meta$term_label, " terms with finite NES and FDR in ", input_file,
    call. = FALSE
  )
}

label_df <- select_labels(plot_df, top_n)
if (nrow(label_df) == 0) {
  message("No significant ", plot_meta$term_label, " terms to label (FDR < ", fdr_cutoff, ").")
}

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

plot_csv <- file.path(RESULTS_DIR, gsea_out_name(paste0(plot_meta$tag, "_volcano_table.csv"), prefix))
write_csv(
  plot_df %>%
    mutate(
      labeled = ID %in% label_df$ID,
      label = as.character(label)
    ),
  plot_csv
)

vol_file <- file.path(RESULTS_DIR, gsea_out_name(paste0(plot_meta$tag, "_volcano.pdf"), prefix))

ggsave(
  vol_file,
  plot = gsea_volcano_plot(plot_df, label_df, plot_meta$title, fdr_cutoff),
  width = 9,
  height = 7
)

n_sig <- sum(plot_df$significant)
dir_counts <- table(label_df$direction)

message("================")
message("Tested ", plot_meta$term_label, " terms: ", nrow(plot_df))
message("Significant ", plot_meta$term_label, " terms (FDR < ", fdr_cutoff, "): ", n_sig)
if (nrow(label_df) > 0) {
  message(
    "Labeled: ",
    paste(paste0(names(dir_counts), "=", as.integer(dir_counts)), collapse = ", ")
  )
}
message("CSV files written:")
message("  ", plot_csv)
message("Plots created:")
message("Volcano plot of ", plot_meta$term_label, " terms: ", vol_file)
