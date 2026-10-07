# Usage:
#   Rscript simplify_GO.R --in gprofiler_GO.csv --outdir res -n 1000 --ont BP --direction none [--organism hsapiens] [--prefix NAME]
#
# Semantic similarity clustering of GO ORA terms with simplifyEnrichment.
# Reads a normalized ORA table (gprofiler, clusterProfiler, fgsea, etc.).
# --prefix is prepended to output file names; default is none.

ora_simplify_options <- function() {
  list(
    usage = "Rscript simplify_GO.R --in FILE --outdir DIR -n 1000 --ont BP --direction none [--organism hsapiens] [--prefix NAME]",
    description = paste(
    "Cluster GO terms by semantic similarity with simplifyEnrichment.",
    "Input is normalized with ora_utils.R so gprofiler, clusterProfiler, and fgsea tables work.",
    "",
    "Examples:",
    "  Rscript simplify_GO.R --in gprofiler_GO.csv --outdir res -n 1000 --ont BP --direction none",
    "  Rscript simplify_GO.R --in gprofiler_GO.csv --outdir res --ont BP --direction up --organism hsapiens",
    "  Rscript simplify_GO.R --in gprofiler_GO.csv --outdir res --ont BP --direction up --prefix young",
    sep = "\n"
  ),
    require_input = TRUE,
    options = list(
    opt_in("gprofiler_GO.csv", "ORA result CSV/TSV or RDS (gprofiler, clusterProfiler, fgsea)"),
    opt_outdir("res", "output directory for PDF plots and tables"),
    opt_n(1000L, "top terms to cluster [default: %default]"),
    opt_ont("BP", metavar = "BP|CC|MF", help = "GO ontology [default: %default]"),
    opt_direction(
      "none",
      metavar = "up|down|none",
      help = paste(
        "term set to cluster: up or down = that direction after filtering,",
        "none = top N terms regardless of direction [default: %default]"
      )
    ),
    opt_organism(
      metavar = "hsapiens|mmusculus",
      help = "organism: hsapiens or mmusculus [default: %default]"
    ),
    opt_prefix()
  )
  )
}

ora_simplify_run <- function(opt) {

  input_file <- opt$input_file
  RESULTS_DIR <- opt$outdir
  top_n <- as.integer(opt$top_n)
  ONT <- toupper(opt$ont)
  direction <- tolower(opt$direction)
  organism <- opt$organism
  prefix <- if (is.null(opt$prefix)) "" else as.character(opt$prefix)
  fdr_cutoff <- 0.05

  organism_to_orgdb <- function(organism) {
    orgdb_map <- c(
      hsapiens = "org.Hs.eg.db",
      mmusculus = "org.Mm.eg.db"
    )
    orgdb <- unname(orgdb_map[organism])
    if (length(orgdb) != 1 || is.na(orgdb)) {
      stop("--organism must be hsapiens or mmusculus", call. = FALSE)
    }
    orgdb
  }

  if (!ONT %in% c("BP", "CC", "MF")) {
    stop("--ont must be BP, CC, or MF", call. = FALSE)
  }
  if (!direction %in% c("up", "down", "none")) {
    stop("--direction must be up, down, or none", call. = FALSE)
  }
  if (!is.finite(top_n) || top_n < 1) {
    stop("-n must be a positive integer", call. = FALSE)
  }
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file, call. = FALSE)
  }

  orgdb_pkg <- organism_to_orgdb(organism)
  if (!requireNamespace(orgdb_pkg, quietly = TRUE)) {
    stop("OrgDb package is not installed: ", orgdb_pkg, call. = FALSE)
  }

  suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(tibble)
    library(simplifyEnrichment)
    library(orgdb_pkg, character.only = TRUE)
  })

  source_sibling("ora_utils.R")

  clusters_to_table <- function(cluster_res, go_ids) {
    if (is.vector(cluster_res) && !is.list(cluster_res) && !is.data.frame(cluster_res)) {
      ids <- names(cluster_res)
      if (is.null(ids) || !any(nzchar(ids))) {
        ids <- go_ids
      }
      return(tibble(term_id = as.character(ids), cluster = as.vector(cluster_res)))
    }

    df <- as.data.frame(cluster_res, stringsAsFactors = FALSE)
    id_col <- intersect(c("id", "term_id", "ID", "go_id", "GO_ID"), colnames(df))
    cl_col <- intersect(c("cluster", "Cluster"), colnames(df))

    if (length(id_col) > 0 && length(cl_col) > 0) {
      return(tibble(
        term_id = as.character(df[[id_col[[1]]]]),
        cluster = df[[cl_col[[1]]]]
      ))
    }

    if (!is.null(rownames(df)) &&
        !identical(rownames(df), as.character(seq_len(nrow(df)))) &&
        length(cl_col) > 0) {
      return(tibble(
        term_id = rownames(df),
        cluster = df[[cl_col[[1]]]]
      ))
    }

    stop(
      "Could not find GO ID and cluster columns in the simplifyGO() result. Columns: ",
      paste(colnames(df), collapse = ", "),
      call. = FALSE
    )
  }

  message("Input:     ", input_file)
  message("Output:    ", RESULTS_DIR)
  message("Prefix:    ", if (nzchar(prefix)) prefix else "(none)")
  message("Organism:  ", organism, " (", orgdb_pkg, ")")
  message("ONT:       ", ONT)
  message("N:         ", top_n)
  message("Direction: ", direction)

  ora_df <- read_ora_plot_table(input_file, ont = ONT, fdr_cutoff = fdr_cutoff)
  if (nrow(ora_df) == 0) {
    stop("No ", ONT, " terms with FDR < ", fdr_cutoff, " in ", input_file, call. = FALSE)
  }

  if (direction %in% c("up", "down")) {
    go_dir <- ora_df %>% filter(.data$direction == .env$direction)
    if (nrow(go_dir) == 0) {
      present <- paste(unique(ora_df$direction), collapse = ", ")
      stop(
        "No ", direction, " ", ONT, " terms after filtering. ",
        "Directions present: ", present, ".",
        call. = FALSE
      )
    }
    go_top <- go_dir %>%
      arrange(p.adjust) %>%
      slice_head(n = top_n)
  } else {
    go_top <- ora_df %>%
      arrange(p.adjust) %>%
      slice_head(n = top_n)
  }

  go_id <- unique(go_top$ID)
  go_id <- go_id[grepl("^GO:[0-9]+$", go_id)]

  if (length(go_id) < 3) {
    stop("Need at least 3 GO IDs for simplifyEnrichment.", call. = FALSE)
  }

  message("GO terms for clustering: ", length(go_id))

  dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

  mat <- GO_similarity(
    go_id,
    ont = ONT,
    db = orgdb_pkg
  )

  stem <- paste0("GO_", ONT, "_", direction, "_simplifyGO")
  plot_file <- file.path(RESULTS_DIR, ora_out_name(paste0(stem, ".pdf"), prefix))
  cluster_csv <- file.path(RESULTS_DIR, ora_out_name(paste0(stem, "_clusters.csv"), prefix))
  sim_rds <- file.path(RESULTS_DIR, ora_out_name(paste0("GO_", ONT, "_", direction, "_similarity.rds"), prefix))

  set.seed(123)
  pdf(plot_file, width = 12, height = 10)
  cluster_res <- simplifyGO(mat)
  dev.off()

  cluster_table <- clusters_to_table(cluster_res, rownames(mat)) %>%
    dplyr::left_join(
      go_top %>%
        transmute(
          term_id = ID,
          term_name = Description,
          p_value = p.adjust,
          intersection_size = Count,
          direction
        ),
      by = "term_id"
    ) %>%
    dplyr::arrange(cluster, p_value)

  write_csv(cluster_table, cluster_csv)
  saveRDS(mat, sim_rds)

  message("================")
  message("Clustered ", nrow(cluster_table), " ", ONT, " terms (direction = ", direction, ")")
  message("CSV files written:")
  message("  ", cluster_csv)
  message("Plots created:")
  message("  ", plot_file)
  message("RDS file written:")
  message("  ", sim_rds)
}

if (sys.nframe() == 0L) {
  source(file.path(dirname(dirname(normalizePath(
    sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)),
    winslash = "/", mustWork = TRUE
  ))), "cli.R"))
  invoke_module(ora_simplify_options, ora_simplify_run)
}
