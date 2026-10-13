# Usage:
#   Rscript clusterProf_KEGG_2.R --in edger.csv --outdir res [--fdr 0.05] [--log2fc 1] [--direction no] [--prefix NAME]
#
# Select significant DEGs by FDR and |log2FC|, map them to Entrez IDs, then
# run clusterProfiler enrichKEGG. With --direction yes, up- and
# down-regulated genes are tested separately; with --direction no
# (default), all significant genes are tested together.
# The background is the measured genes that map to Entrez IDs.
# --prefix is prepended to output file names; default is none.

clusterprof_kegg_options <- function() {
  list(
    usage = "Rscript clusterProf_KEGG_2.R --in FILE --outdir DIR [--fdr 0.05] [--log2fc 1] [--direction no] [--prefix NAME]",
    description = paste(
    "Run KEGG over-representation analysis (ORA) with clusterProfiler::enrichKEGG.",
    "With --direction yes, significant genes are split into up- and down-regulated sets.",
    "With --direction no, all significant genes are tested together.",
    "Genes are mapped from Ensembl to Entrez, then tested against the measured-gene background.",
    "",
    "Examples:",
    "  Rscript clusterProf_KEGG_2.R --in edger.csv --outdir res",
    "  Rscript clusterProf_KEGG_2.R --in edger.csv --outdir res --fdr 0.05 --log2fc 1",
    "  Rscript clusterProf_KEGG_2.R --in edger.csv --outdir res --direction yes",
    "  Rscript clusterProf_KEGG_2.R --in young.csv --outdir res --prefix young",
    sep = "\n"
  ),
    require_input = TRUE,
    options = list(
    opt_in("edger.csv", "edgeR differential expression CSV"),
    opt_outdir("clusterProf_kegg", "output directory for tables and RDS objects"),
    opt_fdr(),
    opt_log2fc(),
    opt_organism(),
    opt_direction(
      "no",
      metavar = "yes/no",
      help = paste(
        "split ORA by regulation direction: yes = up and down separately,",
        "no = all significant genes together [default: %default]"
      )
    ),
    opt_prefix()
  )
  )
}

clusterprof_kegg_run <- function(opt) {

  ORA_FDR <- 0.05

  input_file <- opt$input_file
  RESULTS_DIR <- opt$outdir
  fdr_cutoff <- opt$fdr_cutoff
  log2fc_cutoff <- opt$log2fc_cutoff
  organism <- tolower(opt$organism)
  split_direction <- tolower(opt$direction)
  prefix <- if (is.null(opt$prefix)) "" else as.character(opt$prefix)

  out_name <- function(stem) {
    if (!nzchar(prefix)) {
      return(stem)
    }
    prefix <- gsub("[/\\\\]", "", prefix)
    sep <- if (grepl("[_-]$", prefix)) "" else "_"
    paste0(prefix, sep, stem)
  }

  organism_info <- function(organism) {
    info_map <- list(
      hsapiens = list(orgdb = "org.Hs.eg.db", kegg = "hsa"),
      human = list(orgdb = "org.Hs.eg.db", kegg = "hsa"),
      hs = list(orgdb = "org.Hs.eg.db", kegg = "hsa"),
      mmusculus = list(orgdb = "org.Mm.eg.db", kegg = "mmu"),
      mouse = list(orgdb = "org.Mm.eg.db", kegg = "mmu"),
      mm = list(orgdb = "org.Mm.eg.db", kegg = "mmu")
    )
    info <- info_map[[organism]]
    if (is.null(info)) {
      stop("--organism must be hsapiens or mmusculus", call. = FALSE)
    }
    info
  }

  if (!split_direction %in% c("yes", "no")) {
    stop("--direction must be yes or no", call. = FALSE)
  }
  if (!is.finite(fdr_cutoff) || fdr_cutoff <= 0 || fdr_cutoff > 1) {
    stop("--fdr must be between 0 and 1", call. = FALSE)
  }
  if (!is.finite(log2fc_cutoff) || log2fc_cutoff < 0) {
    stop("--log2fc must be a non-negative number", call. = FALSE)
  }
  if (!file.exists(input_file)) {
    stop("Input file not found: ", input_file, call. = FALSE)
  }

  org <- organism_info(organism)

  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(AnnotationDbi)
    library(readr)
    library(dplyr)
    library(tibble)
  })

  if (!requireNamespace(org$orgdb, quietly = TRUE)) {
    stop("Install the OrgDb package ", org$orgdb, call. = FALSE)
  }
  suppressPackageStartupMessages(library(org$orgdb, character.only = TRUE))
  OrgDb <- get(org$orgdb)

  empty_ora_table <- function() {
    tibble(
      direction = character(),
      ONTOLOGY = character(),
      ID = character(),
      Description = character(),
      p.adjust = double()
    )
  }

  significant_ora_table <- function(df) {
    if (is.null(df) || nrow(df) == 0) {
      return(df)
    }
    df %>% filter(.data$p.adjust < ORA_FDR)
  }

  n_sig_terms <- function(ego) {
    if (is.null(ego)) {
      return(0)
    }
    df <- as.data.frame(ego)
    if (nrow(df) == 0 || !"p.adjust" %in% colnames(df)) {
      return(0)
    }
    sum(df$p.adjust < ORA_FDR, na.rm = TRUE)
  }

  ensembl_to_entrez <- function(ensembl_ids) {
    ensembl_ids <- unique(ensembl_ids)
    ensembl_ids <- ensembl_ids[!is.na(ensembl_ids) & nzchar(ensembl_ids)]
    if (length(ensembl_ids) == 0) {
      return(character())
    }
    mapped <- tryCatch(
      suppressMessages(
        bitr(
          ensembl_ids,
          fromType = "ENSEMBL",
          toType = "ENTREZID",
          OrgDb = OrgDb
        )
      ),
      error = function(e) {
        stop("Ensembl to Entrez mapping failed: ", conditionMessage(e), call. = FALSE)
      }
    )
    unique(mapped$ENTREZID)
  }

  enrich_to_table <- function(ego, direction) {
    if (is.null(ego)) {
      return(NULL)
    }
    df <- as.data.frame(ego)
    if (nrow(df) == 0) {
      return(NULL)
    }
    df$direction <- direction
    df$ONTOLOGY <- "KEGG"
    df %>% relocate(direction, ONTOLOGY, ID, Description, p.adjust)
  }

  run_enrichKEGG <- function(gene_ids, background, direction) {
    gene_ids <- unique(gene_ids)
    gene_ids <- gene_ids[gene_ids %in% background]

    if (length(gene_ids) == 0) {
      message("Skipping ", direction, ": no significant genes mapped to Entrez.")
      return(NULL)
    }

    message(
      "Running clusterProfiler KEGG ORA for ", direction,
      " genes (n = ", length(gene_ids), ")"
    )

    ego <- tryCatch(
      enrichKEGG(
        gene = gene_ids,
        universe = background,
        organism = org$kegg,
        keyType = "kegg",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 10,
        maxGSSize = Inf
      ),
      error = function(e) {
        message("  enrichKEGG failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      return(ego)
    }

    readable <- tryCatch(
      setReadable(ego, OrgDb = OrgDb, keyType = "ENTREZID"),
      error = function(e) {
        message("  setReadable failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(readable)) {
      return(readable)
    }
    ego
  }

  message("Input:     ", input_file)
  message("Output:    ", RESULTS_DIR)
  message("Prefix:    ", if (nzchar(prefix)) prefix else "(none)")
  message("Organism:  ", organism, " (", org$orgdb, ", KEGG ", org$kegg, ")")
  message("DEG FDR:   ", fdr_cutoff)
  message("DEG log2FC:", log2fc_cutoff)
  message("Direction: ", split_direction)

  de <- read_csv(input_file, show_col_types = FALSE)

  required_cols <- c("name", "gene", "log2FoldChange", "FDR")
  missing_cols <- setdiff(required_cols, colnames(de))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # Strip Ensembl version suffixes and keep one row per gene ID.
  gene_map <- de %>%
    mutate(
      ensembl_gene = sub("\\..*$", "", name),
      gene_symbol = as.character(gene)
    ) %>%
    filter(
      !is.na(ensembl_gene),
      ensembl_gene != "",
      !is.na(log2FoldChange),
      is.finite(log2FoldChange),
      !is.na(FDR)
    ) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    distinct(ensembl_gene, .keep_all = TRUE)

  deg_sig <- gene_map %>%
    filter(
      FDR < fdr_cutoff,
      abs(log2FoldChange) >= log2fc_cutoff
    )

  message("Background genes:     ", nrow(gene_map))
  message("Total significant DEGs: ", nrow(deg_sig))

  if (nrow(deg_sig) == 0) {
    stop(
      "No significant genes at FDR < ", fdr_cutoff,
      " and |log2FoldChange| >= ", log2fc_cutoff,
      call. = FALSE
    )
  }

  background <- ensembl_to_entrez(gene_map$ensembl_gene)
  message("Entrez background:    ", length(background))
  if (length(background) == 0) {
    stop("No background genes could be mapped from Ensembl to Entrez.", call. = FALSE)
  }

  dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

  gene_cols <- c("ensembl_gene", "gene_symbol", "log2FoldChange", "PValue", "FDR")
  gene_cols <- gene_cols[gene_cols %in% colnames(gene_map)]

  ora_csv <- file.path(RESULTS_DIR, out_name("clusterProf_KEGG.csv"))
  # ora_all_csv <- file.path(RESULTS_DIR, out_name("clusterProf_KEGG_all.csv"))
  # ora_rds <- file.path(RESULTS_DIR, out_name("clusterProf_KEGG.rds"))

  params <- list(
    input_file = input_file,
    organism = organism,
    orgdb = org$orgdb,
    kegg_organism = org$kegg,
    fdr_cutoff = fdr_cutoff,
    log2fc_cutoff = log2fc_cutoff,
    direction = split_direction,
    prefix = prefix,
    ora_fdr = ORA_FDR,
    keyType = "ENTREZID",
    background_n = length(background)
  )

  gene_csvs <- character()

  if (split_direction == "yes") {
    up_genes <- deg_sig %>%
      filter(log2FoldChange >= log2fc_cutoff)

    down_genes <- deg_sig %>%
      filter(log2FoldChange <= -log2fc_cutoff)

    message("Upregulated:          ", nrow(up_genes))
    message("Downregulated:        ", nrow(down_genes))

    up_genes_csv <- file.path(RESULTS_DIR, out_name("ora_genes_up.csv"))
    down_genes_csv <- file.path(RESULTS_DIR, out_name("ora_genes_down.csv"))
    write_csv(select(up_genes, all_of(gene_cols)), up_genes_csv)
    write_csv(select(down_genes, all_of(gene_cols)), down_genes_csv)
    gene_csvs <- c(up_genes_csv, down_genes_csv)

    enrich_up <- run_enrichKEGG(ensembl_to_entrez(up_genes$ensembl_gene), background, "up")
    enrich_down <- run_enrichKEGG(ensembl_to_entrez(down_genes$ensembl_gene), background, "down")

    result_table <- bind_rows(
      enrich_to_table(enrich_up, "up"),
      enrich_to_table(enrich_down, "down")
    )

    enrich_res <- list(up = enrich_up, down = enrich_down, params = params)
    term_counts <- c(up = n_sig_terms(enrich_up), down = n_sig_terms(enrich_down))
  } else {
    all_genes_csv <- file.path(RESULTS_DIR, out_name("ora_genes_all.csv"))
    write_csv(select(deg_sig, all_of(gene_cols)), all_genes_csv)
    gene_csvs <- all_genes_csv

    enrich_all <- run_enrichKEGG(ensembl_to_entrez(deg_sig$ensembl_gene), background, "none")

    result_table <- enrich_to_table(enrich_all, "none")
    enrich_res <- list(none = enrich_all, params = params)
    term_counts <- c(none = n_sig_terms(enrich_all))
  }

  if (is.null(result_table) || nrow(result_table) == 0) {
    result_table <- empty_ora_table()
  }

  sig_table <- significant_ora_table(result_table)
  if (nrow(sig_table) == 0) {
    if (split_direction == "yes") {
      message("No significant KEGG pathways found for either direction.")
    } else {
      message("No significant KEGG pathways found.")
    }
  }

  write_csv(sig_table, ora_csv)
  # write_csv(result_table, ora_all_csv)
  # saveRDS(enrich_res, ora_rds)

  message("================")
  message("KEGG pathways (clusterProfiler BH FDR < ", ORA_FDR, ")")
  message("  all tested: ", nrow(result_table))
  message("  significant:")
  for (label in names(term_counts)) {
    message("    ", label, ": ", term_counts[[label]])
  }
  message("CSV files written:")
  message("  ", ora_csv)
  # message("  ", ora_all_csv)
  for (gene_csv in gene_csvs) {
    message("  ", gene_csv)
  }
  # message("RDS file written:")
  # message("  ", ora_rds)
}

if (sys.nframe() == 0L) {
  source(file.path(dirname(dirname(normalizePath(
    sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)),
    winslash = "/", mustWork = TRUE
  ))), "cli.R"))
  invoke_module(clusterprof_kegg_options, clusterprof_kegg_run)
}
