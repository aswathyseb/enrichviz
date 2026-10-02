# Usage:
#   Rscript gsea_GO.R --in edger.csv --outdir res [--ont BP] [--GO GO:xyz] [--rank signed_logp]
#
# If --GO is omitted, gseaNb plots the term with the highest positive NES.
# If --rank is omitted, create_gene_list() uses its default (signed_logp).

gsea_go_clusterprof_options <- function() {
  list(
    usage = "Rscript gsea_GO.R --in FILE --outdir DIR [--ont BP] [--GO GO:ID] [--rank METHOD]",
    description = paste(
    "Run GO GSEA with clusterProfiler and save GseaVis plots.",
    "",
    "Examples:",
    "  Rscript gsea_GO.R --in edger.csv --outdir res",
    "  Rscript gsea_GO.R --in edger.csv --outdir res --ont BP --GO GO:0051607 --rank log2FC",
    sep = "\n"
  ),
    require_input = TRUE,
    options = list(
    opt_in("edger.csv", "edgeR differential expression CSV"),
    opt_outdir("gsea_go", "output directory for tables and PDF plots"),
    opt_ont("BP", metavar = "BP|CC|MF", help = "GO ontology [default: %default]"),
    opt_go(metavar = "GO:ID", help = "GO term for gseaNb (default: highest positive NES)"),
    opt_rank("signed_logp", metavar = "signed_logp|log2FC")
  )
  )
}

gsea_go_clusterprof_run <- function(opt) {

  # Gene set size parameters
  MIN_GS_SIZE =100
  MAX_GS_SIZE = 500

  # Database  
  DB = "org.Hs.eg.db"

  # P-value cutoff
  pval_cutoff = 0.05

  # Seed value for reproducibility
  seed_value <- 123

  # Number of terms to plot
  N = 10

  # The input file
  input_file <- opt$input_file
  # The output directory
  RESULTS_DIR <- opt$outdir
  # The ontology type
  ONT <- toupper(opt$ont)
  # The GO term ID
  GO_ID <- opt$GO_ID

  # Optional ranking method; NULL uses create_gene_list() default
  ranking_method <- opt$ranking_method

  if (!ONT %in% c("BP", "CC", "MF")) {
    stop("--ont must be BP, CC, or MF", call. = FALSE)
  }

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
  })

  suppressMessages(suppressWarnings(library(enrichit, quietly = TRUE)))
  suppressMessages(suppressWarnings(library(GseaVis, quietly = TRUE)))

  #
  # Source the shared utility functions.
  # The ranked gene list and expression matrix are created using the utils.R file.
  #

  source_sibling("gsea_utils.R")

  message("Input:  ", input_file)
  message("Output: ", RESULTS_DIR)
  message("ONT:    ", ONT)
  message("Rank:   ", ranking_method)
  if (!is.null(GO_ID)) {
    message("GO_ID:  ", GO_ID)
  }


  #  Read the input file
  de <- read_csv(input_file, show_col_types = FALSE)

  #  Create the ranked gene list
  if (is.null(ranking_method)) {
    ranked <- create_gene_list(de)
  } else {
    ranked <- create_gene_list(de, ranking_method = ranking_method)
  }
  geneList <- ranked$gene_list

  # Parse the expression matrix from the input
  expr <- parse_expr(de)

  # Create the results directory
  dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

  ranked_csv <- file.path(RESULTS_DIR, "ranked_genes_GO.csv")
  write_csv(ranked$rank_df, ranked_csv)
  saveRDS(geneList, file.path(RESULTS_DIR, "ranked_genes_GO.rds"))


  # Set the seed for reproducibility
  set.seed(seed_value)

  # Run the GSEA analysis
  ego <- gseGO(geneList     = geneList,
                OrgDb        = DB,
                keyType      = "ENSEMBL",
                ont          = ONT,
                minGSSize    = MIN_GS_SIZE,
                maxGSSize    = MAX_GS_SIZE,
                pvalueCutoff = pval_cutoff,
                seed         = TRUE,
                verbose      = FALSE)

  # Map Ensembl IDs in the leading-edge / core-enrichment column to gene symbols
  ego <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL")

  # Write the GSEA result table to a file
  result_table <- as.data.frame(ego)
  gsea_csv <- file.path(RESULTS_DIR, paste0("gseGO_", ONT, ".csv"))
  write_csv(result_table, gsea_csv)
  saveRDS(ego, file.path(RESULTS_DIR, paste0("gseGO_", ONT, ".rds")))
  n_sig <- sum(result_table$p.adjust < pval_cutoff, na.rm = TRUE)

  # NES dotplot of top positive and negative terms with dotplotGsea()
  # dot_res <- dotplotGsea(data = ego, topn = N, order.by = "NES", add.seg = TRUE)
  # dot_plot <- if (inherits(dot_res, "ggplot")) dot_res else dot_res$plot

  # ggsave(
  #   file.path(RESULTS_DIR, paste0("gseGO_", ONT, "_dotplot.pdf")),
  #   plot = dot_plot,
  #   width = 10,
  #   height = max(7, 0.3 * N + 2)
  # )

  # NES lollipop of top positive and negative terms
  dot_plot <- nes_lollipop_plot(result_table, ont = ONT, top_n = N, fdr_cutoff = pval_cutoff)
  dot_plot_file <- file.path(RESULTS_DIR, paste0("gseGO_", ONT, "_NES_lollipop.pdf"))

  ggsave(
    dot_plot_file,
    plot = dot_plot,
    width = 9,
    height = max(5, 0.27 * (2 * N) + 2)
  )


  # Create the volcano plot
  vol_res <- volcanoGsea(data = ego)
  vol_plot <- if (inherits(vol_res, "ggplot")) vol_res else vol_res$plot

  vol_plot_file <- file.path(RESULTS_DIR, paste0("gseGO_", ONT, "_volcano.pdf"))
  ggsave(
    vol_plot_file,
    plot = vol_plot,
    width = 8,
    height = 6
  )

  if (is.null(GO_ID)) {
    GO_ID <- top_positive_nes_id(result_table)
  } else if (!GO_ID %in% result_table$ID) {
    stop("GO_ID not found in gseGO results: ", GO_ID, call. = FALSE)
  }

  nb_plot <- gseaNb(object = ego,
                   geneSetID = GO_ID,
                   add.geneExpHt = T,
                   exp = expr)

  nb_file <- file.path(
    RESULTS_DIR,
    paste0("gseGO_", ONT, "_", gsub(":", "_", GO_ID), "_gseaNb.pdf")
  )
  ggsave(
    nb_file,
    plot = nb_plot,
    width = 10,
    height = 10
  )

  message("================")
  message("Significant ", ONT, " terms (p.adjust < ", pval_cutoff, "): ", n_sig)
  message("gseaNb geneSetID: ", GO_ID)

  message("CSV files written:")
  message("  ", ranked_csv)
  message("  ", gsea_csv)
  message("Significant ", ONT, " terms (p.adjust < ", pval_cutoff, "): ", n_sig)
  message("Plots created:")
  message("Lollipop plot of significant ", ONT, " terms: ", dot_plot_file) 
  message("Volcano plot of significant ", ONT, " terms: ", vol_plot_file) 
  message("GSEANES plot of ", GO_ID, ": ", nb_file)
}

if (sys.nframe() == 0L) {
  source(file.path(dirname(dirname(normalizePath(
    sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)),
    winslash = "/", mustWork = TRUE
  ))), "cli.R"))
  invoke_module(gsea_go_clusterprof_options, gsea_go_clusterprof_run)
}
