# Usage:
#   Rscript gsea_KEGG_clusterProf.R --in edger.csv --outdir res [--KEGG hsa04110] [--rank signed_logp]
#
# If --KEGG is omitted, gseaNb plots the pathway with the highest positive NES.
# If --rank is omitted, create_gene_list() uses its default (signed_logp).

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--in", type = "character", default = "edger.csv", dest = "input_file",
              metavar = "FILE", help = "edgeR differential expression CSV"),
  make_option("--outdir", type = "character", default = "gsea_kegg",
              metavar = "DIR", help = "output directory for tables and PDF plots"),
  make_option("--KEGG", type = "character", default = NULL, dest = "KEGG_ID",
              metavar = "ID",
              help = "KEGG pathway for gseaNb (default: highest positive NES)"),
  make_option("--rank", type = "character", default = "signed_logp", dest = "ranking_method",
              metavar = "signed_logp|log2FC",
              help = "ranking statistic [default: signed_logp]")
)

opt_parser <- OptionParser(
  usage = "Rscript gsea_KEGG_clusterProf.R --in FILE --outdir DIR [--KEGG ID] [--rank METHOD]",
  option_list = option_list,
  description = paste(
    "Run KEGG GSEA with clusterProfiler and save GseaVis plots.",
    "",
    "Examples:",
    "  Rscript gsea_KEGG_clusterProf.R --in edger.csv --outdir res",
    "  Rscript gsea_KEGG_clusterProf.R --in edger.csv --outdir res --KEGG hsa04110 --rank log2FC",
    sep = "\n"
  )
)

opt <- parse_args(opt_parser)

if (is.null(opt$input_file) || is.null(opt$outdir)) {
  print_help(opt_parser)
  quit(save = "no", status = 1)
}

MIN_GS_SIZE <- 100
MAX_GS_SIZE <- 500
ORGANISM <- "hsa"
pval_cutoff <- 0.05
seed_value <- 123
N <- 10

input_file <- opt$input_file
RESULTS_DIR <- opt$outdir
KEGG_ID <- opt$KEGG_ID
ranking_method <- opt$ranking_method

if (!is.null(ranking_method) && !ranking_method %in% c("signed_logp", "log2FC")) {
  stop("--rank must be signed_logp or log2FC", call. = FALSE)
}

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

suppressMessages(suppressWarnings(library(enrichit, quietly = TRUE)))
suppressMessages(suppressWarnings(library(GseaVis, quietly = TRUE)))

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

message("Input:  ", input_file)
message("Output: ", RESULTS_DIR)
message("Rank:   ", ranking_method)
if (!is.null(KEGG_ID)) {
  message("KEGG:   ", KEGG_ID)
}

de <- read_csv(input_file, show_col_types = FALSE)

if (is.null(ranking_method)) {
  ranked <- create_gene_list(de)
} else {
  ranked <- create_gene_list(de, ranking_method = ranking_method)
}
geneList_ensembl <- ranked$gene_list
expr <- parse_expr(de)

# gseKEGG() expects Entrez IDs. Map from Ensembl and keep one score per Entrez ID.
id_map <- suppressMessages(
  bitr(
    names(geneList_ensembl),
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
)

rank_entrez <- data.frame(
  ensembl_gene = names(geneList_ensembl),
  rank_score = unname(geneList_ensembl),
  stringsAsFactors = FALSE
) %>%
  inner_join(id_map, by = c("ensembl_gene" = "ENSEMBL")) %>%
  arrange(desc(abs(rank_score))) %>%
  distinct(ENTREZID, .keep_all = TRUE) %>%
  arrange(desc(rank_score))

geneList <- rank_entrez$rank_score
names(geneList) <- rank_entrez$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

if (anyDuplicated(names(geneList)) > 0) {
  stop("Duplicated Entrez IDs remain in the ranked gene list.")
}

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

ranked_csv <- file.path(RESULTS_DIR, "ranked_genes_KEGG.csv")
write_csv(ranked$rank_df, ranked_csv)
saveRDS(geneList, file.path(RESULTS_DIR, "ranked_genes_KEGG.rds"))

set.seed(seed_value)
kk <- gseKEGG(
  geneList = geneList,
  organism = ORGANISM,
  minGSSize = MIN_GS_SIZE,
  maxGSSize = MAX_GS_SIZE,
  pvalueCutoff = pval_cutoff,
  seed = TRUE,
  verbose = FALSE
)

kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

result_table <- as.data.frame(kk)
gsea_csv <- file.path(RESULTS_DIR, "gseKEGG.csv")
write_csv(result_table, gsea_csv)
saveRDS(kk, file.path(RESULTS_DIR, "gseKEGG.rds"))
n_sig <- sum(result_table$p.adjust < pval_cutoff, na.rm = TRUE)

###
dot_res <- dotplotGsea(data = kk, topn = N, order.by = "NES", add.seg = TRUE)
dot_plot <- if (inherits(dot_res, "ggplot")) dot_res else dot_res$plot

ggsave(
  file.path(RESULTS_DIR, paste0("gseKEGG_NES_dotplot.pdf")),
  plot = dot_plot,
  width = 10,
  height = max(7, 0.3 * N + 2)
)
###

# NES lollipop (KEGG-specific; not nes_lollipop_plot from gsea_utils.R)
plot_df <- result_table %>%
  filter(!is.na(p.adjust), p.adjust < pval_cutoff) %>%
  mutate(
    direction = ifelse(NES > 0, "up", "down"),
    label = str_to_sentence(Description)
  ) %>%
  group_by(direction) %>%
  slice_max(order_by = abs(NES), n = N, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(NES) %>%
  mutate(label = factor(label, levels = label))

if (nrow(plot_df) == 0) {
  stop("No significant KEGG pathways to plot in the NES lollipop.")
}

dot_plot <- ggplot(plot_df, aes(x = NES, y = label, color = p.adjust)) +
  geom_segment(aes(x = 0, xend = NES, yend = label), linewidth = 0.6, color = "grey70") +
  geom_point(aes(size = setSize), alpha = 0.9) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
  scale_color_viridis_c(option = "magma", direction = -1, trans = "log10") +
  labs(
    title = "KEGG GSEA",
    subtitle = "Top positively and negatively enriched pathways",
    x = "Normalized enrichment score (NES)",
    y = NULL,
    color = "FDR",
    size = "Gene-set size"
  ) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 10))

dot_plot_file <- file.path(RESULTS_DIR, "gseKEGG_NES_lollipop.pdf")
ggsave(
  dot_plot_file,
  plot = dot_plot,
  width = 9,
  height = max(5, 0.27 * nrow(plot_df) + 2)
)

vol_res <- volcanoGsea(data = kk)
vol_plot <- if (inherits(vol_res, "ggplot")) vol_res else vol_res$plot
vol_plot_file <- file.path(RESULTS_DIR, "gseKEGG_volcano.pdf")
ggsave(
  vol_plot_file,
  plot = vol_plot,
  width = 8,
  height = 6
)

if (is.null(KEGG_ID)) {
  KEGG_ID <- top_positive_nes_id(result_table)
} else if (!KEGG_ID %in% result_table$ID) {
  stop("KEGG ID not found in gseKEGG results: ", KEGG_ID, call. = FALSE)
}

nb_plot <- gseaNb(
  object = kk,
  geneSetID = KEGG_ID,
  add.geneExpHt = TRUE,
  exp = expr,
  kegg = TRUE
)

nb_file <- file.path(
  RESULTS_DIR,
  paste0("gseKEGG_", gsub("[^A-Za-z0-9]+", "_", KEGG_ID), "_gseaNb.pdf")
)
ggsave(
  nb_file,
  plot = nb_plot,
  width = 10,
  height = 10
)

message("================")
message("Significant KEGG pathways (p.adjust < ", pval_cutoff, "): ", n_sig)
message("gseaNb geneSetID: ", KEGG_ID)
message("CSV files written:")
message("  ", ranked_csv)
message("  ", gsea_csv)
message("Plots created:")
message("Lollipop plot of significant KEGG pathways: ", dot_plot_file)
message("Volcano plot of significant KEGG pathways: ", vol_plot_file)
message("GSEANES plot of ", KEGG_ID, ": ", nb_file)
