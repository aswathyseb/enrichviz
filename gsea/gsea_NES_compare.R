# Usage:
#   Rscript gsea_NES_compare.R --GO GO:0051607 \
#     --in young=gsea_out/fgsea_GO_all.rds,old=gsea_out1/gseGO_BP.rds \
#     --outdir gsea_out [--prefix NAME]
#   Rscript gsea_NES_compare.R --term hsa03010 \
#     --in a=t5/fgsea_KEGG.rds,b=t5/fgsea_KEGG.rds --outdir gsea_out
#
# Overlay running-enrichment curves for one GO or KEGG term across conditions.
# Each --in file should be an fgsea/clusterProfiler RDS (preferred) or a
# result table plus a matching --rank file so the ranked list is available.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--in", type = "character", default = NULL, dest = "input_spec",
              metavar = "FILES",
              help = "comma-separated GSEA RDS/CSV files, optionally name=file"),
  make_option("--labels", type = "character", default = NULL,
              metavar = "NAMES",
              help = "comma-separated condition names (if --in files are unlabeled)"),
  make_option("--rank", type = "character", default = NULL, dest = "rank_spec",
              metavar = "FILES",
              help = "optional ranked-gene files, same order or name=file"),
  make_option("--GO", type = "character", default = NULL, dest = "GO_ID",
              metavar = "ID", help = "GO or KEGG term to plot (required unless --term)"),
  make_option("--term", type = "character", default = NULL, dest = "term_id",
              metavar = "ID", help = "alias for --GO"),
  make_option("--outdir", type = "character", default = "gsea_out",
              metavar = "DIR", help = "output directory for PDF plots"),
  make_option("--prefix", type = "character", default = "",
              metavar = "NAME",
              help = "optional prefix for output file names [default: none]")
)

opt_parser <- OptionParser(
  usage = "\nRscript gsea_NES_compare.R --GO ID --in name=FILE,... --outdir DIR [--prefix NAME]\n",
  option_list = option_list,
  description = paste(
    "Compare running enrichment of one GO or KEGG term across GSEA result files.",
    "Each curve is that condition's ranked list walked for the same term.",
    "",
    "Examples:",
    "  Rscript gsea_NES_compare.R --GO GO:0051607 --in young=gsea_out/fgsea_GO_all.rds,old=gsea_out1/gseGO_BP.rds --outdir gsea_out",
    "  Rscript gsea_NES_compare.R --term hsa03010 --in a=t5/fgsea_KEGG.rds,b=t5/fgsea_KEGG.rds --outdir gsea_out --prefix kegg",
    sep = "\n"
  )
)

opt <- parse_args(opt_parser)

term_id <- opt$term_id
if (is.null(term_id) || !nzchar(term_id)) {
  term_id <- opt$GO_ID
}

if (is.null(opt$input_spec) || is.null(opt$outdir) || is.null(term_id) || !nzchar(term_id)) {
  print_help(opt_parser)
  quit(save = "no", status = 1)
}

GO_ID <- term_id
RESULTS_DIR <- opt$outdir
prefix <- if (is.null(opt$prefix)) "" else as.character(opt$prefix)

parse_named_files <- function(spec, labels_spec = NULL, arg_name = "--in") {
  parts <- trimws(unlist(strsplit(spec, ",", fixed = TRUE)))
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) {
    stop(arg_name, " did not contain any files.", call. = FALSE)
  }

  has_eq <- grepl("=", parts, fixed = TRUE)
  if (all(has_eq)) {
    labs <- trimws(sub("=.*$", "", parts))
    files <- trimws(sub("^[^=]*=", "", parts))
  } else if (!any(has_eq)) {
    files <- parts
    if (!is.null(labels_spec) && nzchar(labels_spec)) {
      labs <- trimws(unlist(strsplit(labels_spec, ",", fixed = TRUE)))
      if (length(labs) != length(files)) {
        stop("--labels must have one name per file in ", arg_name, ".", call. = FALSE)
      }
    } else {
      labs <- tools::file_path_sans_ext(basename(files))
    }
  } else {
    stop("Use either name=file,... or a plain file list with --labels.", call. = FALSE)
  }

  if (anyDuplicated(labs) > 0) {
    stop("Condition names must be unique.", call. = FALSE)
  }

  tibble::tibble(condition = labs, file = files)
}

split_rank_files <- function(rank_spec, conditions) {
  ranks <- setNames(rep(NA_character_, length(conditions)), conditions)
  if (is.null(rank_spec) || !nzchar(rank_spec)) {
    return(ranks)
  }

  spec <- parse_named_files(rank_spec, arg_name = "--rank")
  if (all(spec$condition %in% conditions)) {
    ranks[spec$condition] <- spec$file
  } else if (nrow(spec) == length(conditions)) {
    ranks[] <- spec$file
  } else {
    stop("--rank must be name=file for each condition, or a list in the same order as --in.", call. = FALSE)
  }
  ranks
}

cond_tbl <- parse_named_files(opt$input_spec, opt$labels)
if (nrow(cond_tbl) < 2) {
  stop("Provide at least two GSEA result files in --in to compare.", call. = FALSE)
}
rank_map <- split_rank_files(opt$rank_spec, cond_tbl$condition)

missing <- cond_tbl$file[!file.exists(cond_tbl$file)]
if (length(missing) > 0) {
  stop("Input file not found: ", paste(missing, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
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

COMPARE_COLORS <- c(
  "#4DAF4A", "#2C4A7C", "#E07A3D", "#984EA3",
  "#E41A1C", "#377EB8", "#A65628", "#F781BF"
)

load_condition <- function(label, path, rank_file, go_id) {
  rank_arg <- if (!is.na(rank_file) && nzchar(rank_file)) rank_file else NULL
  nes_in <- tryCatch(
    load_gsea_nes_term(path, go_id, rank_file = rank_arg),
    error = function(e) {
      message("Skipping ", label, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(nes_in) || is.null(nes_in$stats)) {
    return(NULL)
  }
  es_df <- gsea_running_es(nes_in$stats, nes_in$pathway_genes)
  es_df$condition <- label
  list(es = es_df, term = nes_in$term, go_id = nes_in$go_id)
}

message("Term:   ", GO_ID)
message("Output: ", RESULTS_DIR)
for (i in seq_len(nrow(cond_tbl))) {
  message("  ", cond_tbl$condition[[i]], ": ", cond_tbl$file[[i]])
}

loaded <- Map(
  load_condition,
  cond_tbl$condition,
  cond_tbl$file,
  rank_map[cond_tbl$condition],
  MoreArgs = list(go_id = GO_ID)
)
loaded <- loaded[!vapply(loaded, is.null, logical(1))]
if (length(loaded) < 2) {
  stop("Need at least two conditions with a ranked list and the selected term.", call. = FALSE)
}

es_all <- bind_rows(lapply(loaded, `[[`, "es"))
cond_levels <- cond_tbl$condition[cond_tbl$condition %in% unique(es_all$condition)]
es_all$condition <- factor(es_all$condition, levels = cond_levels)

summary_tbl <- bind_rows(lapply(names(loaded), function(nm) {
  term <- loaded[[nm]]$term
  tibble::tibble(
    condition = nm,
    ID = term$ID[[1]],
    Description = term$Description[[1]],
    NES = term$NES[[1]],
    p.adjust = term$p.adjust[[1]],
    n_in_set = sum(loaded[[nm]]$es$in_set),
    n_ranked = nrow(loaded[[nm]]$es)
  )
}))

title <- str_to_title(summary_tbl$Description[[1]])
xmax <- max(es_all$rank)
n_cond <- nlevels(es_all$condition)
pal <- setNames(COMPARE_COLORS[seq_len(n_cond)], cond_levels)

ticks <- es_all %>%
  filter(in_set) %>%
  mutate(condition = factor(condition, levels = rev(cond_levels)))

es_plot <- ggplot(es_all, aes(x = rank, y = running_es, color = condition, fill = condition)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35, color = "grey50") +
  geom_ribbon(
    aes(ymin = pmin(running_es, 0), ymax = pmax(running_es, 0)),
    alpha = 0.18,
    color = NA
  ) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(sec.axis = dup_axis(name = "ES")) +
  coord_cartesian(xlim = c(0, xmax), expand = FALSE) +
  labs(
    title = title,
    x = "Rank",
    y = "ES",
    color = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "inside",
    legend.position.inside = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.4, "cm")
  )

barcode_plot <- ggplot(ticks, aes(x = rank, y = condition, color = condition)) +
  geom_segment(
    aes(xend = rank, y = as.numeric(condition) - 0.38, yend = as.numeric(condition) + 0.38),
    linewidth = 0.25,
    alpha = 0.9
  ) +
  scale_color_manual(values = pal, guide = "none") +
  scale_y_discrete(drop = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, xmax), expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

grad_df <- tibble::tibble(
  x = seq(0, xmax, length.out = 400),
  y = 1
)
grad_plot <- ggplot(grad_df, aes(x = x, y = y, fill = x)) +
  geom_raster() +
  scale_fill_gradientn(colors = c("#B2182B", "#F7F7F7", "#2166AC"), guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, xmax), expand = FALSE) +
  annotate(
    "text",
    x = xmax * 0.02,
    y = 1,
    label = "positively correlated",
    hjust = 0,
    vjust = 0.5,
    size = 3,
    color = "grey20"
  ) +
  annotate(
    "text",
    x = xmax * 0.98,
    y = 1,
    label = "negatively correlated",
    hjust = 1,
    vjust = 0.5,
    size = 3,
    color = "grey20"
  ) +
  labs(x = NULL, y = NULL) +
  theme_void()

barcode_h <- max(0.55, 0.28 * n_cond)
combined <- es_plot / barcode_plot / grad_plot +
  plot_layout(heights = c(3.4, barcode_h, 0.32))

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
if (nzchar(prefix)) message("Prefix: ", prefix)

go_tag <- gsea_term_tag(GO_ID)
plot_csv <- file.path(RESULTS_DIR, gsea_out_name(paste0("gsea_", go_tag, "_NES_compare_table.csv"), prefix))
write_csv(summary_tbl, plot_csv)

plot_file <- file.path(RESULTS_DIR, gsea_out_name(paste0("gsea_", go_tag, "_NES_compare.pdf"), prefix))
ggsave(plot_file, plot = combined, width = 8.5, height = 5.4 + barcode_h)

message("================")
message("Term: ", summary_tbl$ID[[1]], " — ", summary_tbl$Description[[1]])
for (i in seq_len(nrow(summary_tbl))) {
  message(
    "  ", summary_tbl$condition[[i]],
    ": NES=", signif(summary_tbl$NES[[i]], 4),
    "  FDR=", signif(summary_tbl$p.adjust[[i]], 3)
  )
}
message("CSV files written:")
message("  ", plot_csv)
message("Plots created:")
message("NES compare plot: ", plot_file)
