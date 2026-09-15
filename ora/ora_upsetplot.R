# Usage:
#   Rscript ora_upsetplot.R --in gprofiler_GO.csv --outdir res -n 10 --ont BP --direction up
#   Rscript ora_upsetplot.R --in kegg.csv --outdir res -n 10 --direction none [--prefix NAME]
#
# Upset plot of gene overlap among top ORA terms from gprofiler,
# clusterProfiler, fgsea, or a similar table (GO, KEGG, or another
# gene-set database). --direction selects which terms to plot: up, down,
# or none (top N overall).

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--in", type = "character", default = "gprofiler_GO.csv", dest = "input_file",
              metavar = "FILE",
              help = "ORA result CSV/TSV or RDS (gprofiler, clusterProfiler, fgsea)"),
  make_option("--outdir", type = "character", default = "res",
              metavar = "DIR", help = "output directory for PDF plots"),
  make_option(c("-n", "--n"), type = "integer", default = 10, dest = "top_n",
              metavar = "N", help = "top terms to plot [default: %default]"),
  make_option("--ont", type = "character", default = "BP",
              metavar = "BP|CC|MF|KEGG|all",
              help = paste(
                "GO ontology to plot, or KEGG/all to skip GO filtering",
                "[default: %default]"
              )),
  make_option("--direction", type = "character", default = "none",
              metavar = "up|down|none",
              help = "terms to plot: up, down, or none (top N overall) [default: %default]"),
  make_option("--prefix", type = "character", default = "",
              metavar = "NAME",
              help = "optional prefix for output file names [default: none]")
)

opt_parser <- OptionParser(
  usage = "\nRscript ora_upsetplot.R --in FILE --outdir DIR -n 10 [--ont BP] [--direction up] [--prefix NAME]",
  option_list = option_list,
  description = paste(
    "\nUpset plot of gene overlap among top ORA terms from GO, KEGG, or another gene-set database.\n",
    "Accepts gprofiler, clusterProfiler, or fgsea tables (CSV/TSV/RDS).",
    "For GO tables, --ont BP, CC, or MF selects one ontology.",
    "For KEGG and other non-GO tables, ontology filtering is skipped. \n",
    "Use --direction up or down to plot that subset; none plots the top N terms overall.",
    "--prefix is prepended to output file names.",
    "",
    "Examples:",
    "  Rscript ora_upsetplot.R --in gprofiler_GO.csv --outdir res -n 10 --ont BP --direction up",
    "  Rscript ora_upsetplot.R --in kegg.csv --outdir res -n 10 --direction none --prefix KEGG",
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
DIRECTION <- tolower(opt$direction)
prefix <- if (is.null(opt$prefix)) "" else as.character(opt$prefix)
fdr_cutoff <- 0.05

if (!ONT %in% c("BP", "CC", "MF", "KEGG", "ALL")) {
  stop("--ont must be BP, CC, MF, KEGG, or all", call. = FALSE)
}
if (!DIRECTION %in% c("up", "down", "none")) {
  stop("--direction must be up, down, or none", call. = FALSE)
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
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(ComplexUpset)
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

select_top_terms <- function(df, n, direction) {
  if (!"direction" %in% colnames(df) || all(is.na(df$direction) | df$direction == "")) {
    df$direction <- "none"
  }

  if (direction %in% c("up", "down")) {
    df <- df %>% filter(.data$direction == !!direction)
  }

  plot_df <- df %>%
    slice_min(order_by = p.adjust, n = n, with_ties = FALSE) %>%
    mutate(label = str_to_sentence(Description)) %>%
    arrange(p.adjust)

  if (anyDuplicated(plot_df$label) > 0) {
    plot_df <- plot_df %>%
      mutate(label = paste0(label, " (", ID, ")"))
  }

  plot_df %>%
    mutate(label = factor(label, levels = unique(label)))
}

build_gene_membership <- function(plot_df) {
  membership <- plot_df %>%
    transmute(label = as.character(label), geneID = as.character(geneID)) %>%
    filter(!is.na(geneID), geneID != "") %>%
    separate_rows(geneID, sep = "[,/;]") %>%
    mutate(geneID = str_trim(geneID)) %>%
    filter(geneID != "") %>%
    distinct(label, geneID) %>%
    mutate(value = TRUE) %>%
    pivot_wider(
      names_from = label,
      values_from = value,
      values_fill = FALSE
    )

  membership
}

ora_upset_plot <- function(membership, set_names) {
  upset(
    membership,
    intersect = set_names,
    set_sizes = FALSE
  )
}

message("Input:  ", input_file)
message("Output: ", RESULTS_DIR)
message("Prefix: ", if (nzchar(prefix)) prefix else "(none)")
message("N:      ", top_n)
message("Dir:    ", DIRECTION)

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

plot_df <- select_top_terms(ora_df, top_n, DIRECTION)
if (nrow(plot_df) == 0) {
  if (DIRECTION %in% c("up", "down")) {
    stop(
      "No ", DIRECTION, "-regulated ", plot_meta$term_label, " terms with FDR < ", fdr_cutoff, ".",
      call. = FALSE
    )
  }
  stop("No terms left to plot after selecting the top ", top_n, ".", call. = FALSE)
}

membership <- build_gene_membership(plot_df)
set_names <- intersect(as.character(plot_df$label), colnames(membership))
if (nrow(membership) == 0 || length(set_names) < 2) {
  stop(
    "Need gene lists for at least two terms to draw an upset plot.",
    call. = FALSE
  )
}

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

if (DIRECTION == "none") {
  stem <- paste0(plot_meta$tag, "_upset")
} else {
  stem <- paste0(plot_meta$tag, "_", DIRECTION, "_upset")
}
plot_csv <- file.path(RESULTS_DIR, ora_out_name(paste0(stem, "_table.csv"), prefix))
write_csv(plot_df %>% mutate(label = as.character(label)), plot_csv)

n_sets <- length(set_names)
plot_width <- max(10, min(20, 0.8 * n_sets + 6))
plot_height <- max(6, min(14, 0.45 * n_sets + 4))

upset_file <- file.path(RESULTS_DIR, ora_out_name(paste0(stem, ".pdf"), prefix))

ggsave(
  upset_file,
  plot = ora_upset_plot(membership, set_names),
  width = plot_width,
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
message("Upset plot of significant ", plot_meta$term_label, " terms: ", upset_file)
