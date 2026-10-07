# Usage:
#   Rscript ora_enricher.R --in edger.csv --outdir res [--fdr 0.05] [--log2fc 1] [--direction no] [--gmt FILE] [--prefix NAME]
#
# Select significant DEGs by FDR and |log2FC|, then run clusterProfiler::enricher.
# This is the local hypergeometric ORA used by gseapy.enrich.
# With --direction yes, up- and down-regulated genes are tested separately.
# With --direction no (default), all significant genes are tested together.
# Without --gmt, gene sets are MSigDB C5 GO for BP, MF, and CC.
# --gmt supplies a custom GMT library instead.
# The background is the set of measured genes in the DE table.
# --prefix is prepended to output file names; default is none.

ora_enricher_options <- function() {
  list(
    usage = "Rscript ora_enricher.R --in FILE --outdir DIR [--fdr 0.05] [--log2fc 1] [--direction no] [--gmt FILE] [--prefix NAME]",
    description = paste(
    "Run over-representation analysis (ORA) with clusterProfiler::enricher.",
    "This is the R counterpart of gseapy.enrich: a local hypergeometric test",
    "against a gene-set library, with the measured genes as the background.",
    "With --direction yes, significant genes are split into up- and down-regulated sets.",
    "With --direction no, all significant genes are tested together.",
    "Without --gmt, MSigDB GO gene sets for BP, MF, and CC are used.",
    "With --gmt, that GMT file is the gene-set library.",
    "",
    "Examples:",
    "  Rscript ora_enricher.R --in edger.csv --outdir res",
    "  Rscript ora_enricher.R --in edger.csv --outdir res --fdr 0.05 --log2fc 1",
    "  Rscript ora_enricher.R --in edger.csv --outdir res --direction yes",
    "  Rscript ora_enricher.R --in edger.csv --outdir res --gmt genesets.gmt --prefix custom",
    sep = "\n"
  ),
    require_input = TRUE,
    options = list(
    opt_in("edger.csv", "edgeR differential expression CSV"),
    opt_outdir("enricher", "output directory for tables and RDS objects"),
    opt_fdr(),
    opt_log2fc(),
    opt_organism(help = "organism for MSigDB GO gene sets: hsapiens or mmusculus [default: %default]"),
    make_option(
      "--gmt",
      type = "character",
      default = NULL,
      dest = "gmt",
      metavar = "FILE",
      help = "GMT gene-set file. When omitted, MSigDB GO gene sets for BP, CC, and MF are used"
    ),
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

ora_enricher_run <- function(opt) {

  GO_ONTS <- c("BP", "MF", "CC")
  ORA_FDR <- 0.05
  MIN_GS_SIZE <- 10
  MAX_GS_SIZE <- 500

  input_file <- opt$input_file
  RESULTS_DIR <- opt$outdir
  fdr_cutoff <- opt$fdr_cutoff
  log2fc_cutoff <- opt$log2fc_cutoff
  organism <- tolower(opt$organism)
  gmt_file <- opt$gmt
  split_direction <- tolower(opt$direction)
  prefix <- if (is.null(opt$prefix)) "" else as.character(opt$prefix)

  if (!is.null(gmt_file) && !nzchar(gmt_file)) {
    gmt_file <- NULL
  }
  use_gmt <- !is.null(gmt_file)

  out_name <- function(stem) {
    if (!nzchar(prefix)) {
      return(stem)
    }
    prefix <- gsub("[/\\\\]", "", prefix)
    sep <- if (grepl("[_-]$", prefix)) "" else "_"
    paste0(prefix, sep, stem)
  }

  organism_to_msig <- function(organism) {
    org_map <- list(
      hsapiens = list(species = "Homo sapiens", db_species = "HS"),
      human = list(species = "Homo sapiens", db_species = "HS"),
      hs = list(species = "Homo sapiens", db_species = "HS"),
      mmusculus = list(species = "Mus musculus", db_species = "MM"),
      mouse = list(species = "Mus musculus", db_species = "MM"),
      mm = list(species = "Mus musculus", db_species = "MM")
    )
    info <- org_map[[organism]]
    if (is.null(info)) {
      stop("--organism must be hsapiens or mmusculus", call. = FALSE)
    }
    info
  }

  go_term_name <- function(gs_name) {
    stripped <- sub("^GO(BP|CC|MF)_", "", gs_name)
    str_to_sentence(gsub("_", " ", stripped))
  }

  read_gmt_sets <- function(path) {
    lines <- readLines(path, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    if (length(lines) == 0) {
      stop("GMT file is empty: ", path, call. = FALSE)
    }
    rows <- strsplit(lines, "\t", fixed = TRUE)

    gene_rows <- lapply(rows, function(x) {
      if (length(x) < 3) {
        return(NULL)
      }
      genes <- x[-c(1, 2)]
      genes <- genes[nzchar(genes)]
      if (length(genes) == 0 || !nzchar(x[[1]])) {
        return(NULL)
      }
      data.frame(term = x[[1]], gene = genes, stringsAsFactors = FALSE)
    })
    term2gene <- bind_rows(gene_rows)
    if (nrow(term2gene) == 0) {
      stop(
        "No gene sets found in ", path,
        ". Expected GMT rows: name, description, gene, gene, ...",
        call. = FALSE
      )
    }

    term2name <- bind_rows(lapply(rows, function(x) {
      if (length(x) < 1 || !nzchar(x[[1]])) {
        return(NULL)
      }
      desc <- if (length(x) >= 2) x[[2]] else ""
      if (!nzchar(desc) || grepl("^https?://", desc)) {
        desc <- str_to_sentence(gsub("_", " ", x[[1]]))
      }
      data.frame(term = x[[1]], name = desc, stringsAsFactors = FALSE)
    })) %>%
      distinct(.data$term, .keep_all = TRUE)

    list(
      term2gene = distinct(term2gene, .data$term, .data$gene),
      term2name = term2name,
      ont_map = setNames(rep(NA_character_, nrow(term2name)), term2name$term),
      id_type = NULL,
      label = basename(path)
    )
  }

  load_go_sets <- function(onts, msig_org) {
    pieces <- lapply(onts, function(ont) {
      subcollection <- paste0("GO:", ont)
      message("Loading MSigDB ", subcollection, " gene sets...")
      msig <- tryCatch(
        msigdbr(
          db_species = msig_org$db_species,
          species = msig_org$species,
          collection = "C5",
          subcollection = subcollection
        ),
        error = function(e) {
          stop("Failed to load MSigDB gene sets: ", conditionMessage(e), call. = FALSE)
        }
      )
      msig <- msig %>%
        filter(
          !is.na(ensembl_gene), ensembl_gene != "",
          !is.na(gs_exact_source), gs_exact_source != "",
          !is.na(gs_name), gs_name != ""
        )
      if (nrow(msig) == 0) {
        stop("No MSigDB ", subcollection, " gene sets with Ensembl IDs were returned.", call. = FALSE)
      }
      message("  gene sets: ", dplyr::n_distinct(msig$gs_exact_source))
      list(
        term2gene = msig %>%
          distinct(.data$gs_exact_source, .data$ensembl_gene) %>%
          transmute(term = .data$gs_exact_source, gene = .data$ensembl_gene),
        term2name = msig %>%
          distinct(.data$gs_exact_source, .data$gs_name) %>%
          transmute(term = .data$gs_exact_source, name = go_term_name(.data$gs_name))
      )
    })

    term2gene <- bind_rows(lapply(pieces, `[[`, "term2gene")) %>%
      distinct(.data$term, .data$gene)
    term2name <- bind_rows(lapply(pieces, `[[`, "term2name")) %>%
      distinct(.data$term, .keep_all = TRUE)
    ont_map <- unlist(lapply(seq_along(pieces), function(i) {
      setNames(
        rep(onts[[i]], nrow(pieces[[i]]$term2name)),
        pieces[[i]]$term2name$term
      )
    }), use.names = TRUE)

    list(
      term2gene = term2gene,
      term2name = term2name,
      ont_map = ont_map,
      id_type = "ensembl",
      label = paste("MSigDB C5", paste(onts, collapse = "+"))
    )
  }

  choose_id_type <- function(gmt_genes, gene_map) {
    gmt_genes <- unique(gmt_genes)
    n_ens <- sum(gmt_genes %in% gene_map$ensembl_gene)
    n_sym <- sum(toupper(gmt_genes) %in% toupper(gene_map$gene_symbol))
    if (n_ens == 0 && n_sym == 0) {
      stop(
        "GMT genes match neither Ensembl IDs nor gene symbols in the DE table.",
        call. = FALSE
      )
    }
    if (n_sym > n_ens) "symbol" else "ensembl"
  }

  ids_to_symbols <- function(gene_id, gene_map) {
    vapply(strsplit(as.character(gene_id), "/", fixed = TRUE), function(ids) {
      ids <- ids[nzchar(ids)]
      if (length(ids) == 0) {
        return("")
      }
      syms <- gene_map$gene_symbol[match(ids, gene_map$ensembl_gene)]
      paste(ifelse(is.na(syms) | syms == "", ids, syms), collapse = "/")
    }, character(1))
  }

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
    df %>% filter(!is.na(.data$p.adjust), .data$p.adjust < ORA_FDR)
  }

  n_sig_terms <- function(df) {
    if (is.null(df) || nrow(df) == 0 || !"p.adjust" %in% colnames(df)) {
      return(0)
    }
    sum(df$p.adjust < ORA_FDR, na.rm = TRUE)
  }

  enrich_to_table <- function(ego, direction, ont_map, id_type, gene_map) {
    if (is.null(ego)) {
      return(NULL)
    }
    df <- as.data.frame(ego)
    if (nrow(df) == 0) {
      return(NULL)
    }
    df$direction <- direction
    df$ONTOLOGY <- unname(ont_map[as.character(df$ID)])
    if (identical(id_type, "ensembl") && "geneID" %in% colnames(df)) {
      df$geneID <- ids_to_symbols(df$geneID, gene_map)
    }
    df %>% relocate(direction, ONTOLOGY, ID, Description, p.adjust)
  }

  run_enricher <- function(gene_ids, background, sets, direction) {
    gene_ids <- unique(gene_ids)
    gene_ids <- gene_ids[gene_ids %in% background]
    if (length(gene_ids) == 0) {
      message("Skipping ", direction, ": no significant genes at the chosen cutoffs.")
      return(NULL)
    }

    if (direction == "none") {
      message("Running enricher ORA for all genes (n = ", length(gene_ids), ")")
    } else {
      message("Running enricher ORA for ", direction, " genes (n = ", length(gene_ids), ")")
    }

    # enricher indexes TERM2GENE with [,], which must drop to an atomic vector.
    term2gene <- as.data.frame(sets$term2gene, stringsAsFactors = FALSE)
    term2name <- as.data.frame(sets$term2name, stringsAsFactors = FALSE)

    tryCatch(
      enricher(
        gene = gene_ids,
        universe = background,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = MIN_GS_SIZE,
        maxGSSize = MAX_GS_SIZE
      ),
      error = function(e) {
        message("  enricher failed: ", conditionMessage(e))
        NULL
      }
    )
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
  if (use_gmt && !file.exists(gmt_file)) {
    stop("GMT file not found: ", gmt_file, call. = FALSE)
  }

  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(readr)
    library(dplyr)
    library(tibble)
    library(stringr)
  })
  if (!use_gmt) {
    suppressPackageStartupMessages(library(msigdbr))
  }

  message("Input:     ", input_file)
  message("Output:    ", RESULTS_DIR)
  message("Prefix:    ", if (nzchar(prefix)) prefix else "(none)")
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

  if (nrow(deg_sig) == 0) {
    stop(
      "No significant genes at FDR < ", fdr_cutoff,
      " and |log2FoldChange| >= ", log2fc_cutoff,
      call. = FALSE
    )
  }

  if (use_gmt) {
    sets <- read_gmt_sets(gmt_file)
    sets$id_type <- choose_id_type(sets$term2gene$gene, gene_map)
    if (sets$id_type == "symbol") {
      sets$term2gene$gene <- toupper(sets$term2gene$gene)
      sets$term2gene <- distinct(sets$term2gene, .data$term, .data$gene)
    }
    message("Gene sets: ", sets$label, " (", sets$id_type, " ids)")
  } else {
    msig_org <- organism_to_msig(organism)
    message(
      "Organism:  ", organism,
      " (", msig_org$species, ", MSigDB ", msig_org$db_species, ")"
    )
    message("Ontology:  ", paste(GO_ONTS, collapse = ", "))
    sets <- load_go_sets(GO_ONTS, msig_org)
  }

  if (sets$id_type == "ensembl") {
    background <- unique(gene_map$ensembl_gene)
    query_ids <- function(df) df$ensembl_gene
  } else {
    background <- unique(toupper(gene_map$gene_symbol))
    background <- background[!is.na(background) & nzchar(background)]
    query_ids <- function(df) {
      ids <- toupper(df$gene_symbol)
      ids[!is.na(ids) & nzchar(ids)]
    }
  }

  message("Background genes:       ", length(background))
  message("Total significant DEGs: ", nrow(deg_sig))
  message("Library size:           ", dplyr::n_distinct(sets$term2gene$term), " gene sets")

  dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

  gene_cols <- c("ensembl_gene", "gene_symbol", "log2FoldChange", "PValue", "FDR")
  gene_cols <- gene_cols[gene_cols %in% colnames(gene_map)]

  stem <- if (use_gmt) "enricher" else "enricher_GO"
  ora_csv <- file.path(RESULTS_DIR, out_name(paste0(stem, ".csv")))
  ora_all_csv <- file.path(RESULTS_DIR, out_name(paste0(stem, "_all.csv")))
  ora_rds <- file.path(RESULTS_DIR, out_name(paste0(stem, ".rds")))

  params <- list(
    input_file = input_file,
    organism = organism,
    fdr_cutoff = fdr_cutoff,
    log2fc_cutoff = log2fc_cutoff,
    direction = split_direction,
    prefix = prefix,
    ora_fdr = ORA_FDR,
    gmt = gmt_file,
    gene_sets = sets$label,
    id_type = sets$id_type,
    minGSSize = MIN_GS_SIZE,
    maxGSSize = MAX_GS_SIZE,
    background_n = length(background)
  )

  gene_csvs <- character()

  if (split_direction == "yes") {
    up_genes <- deg_sig %>%
      filter(log2FoldChange >= log2fc_cutoff)

    down_genes <- deg_sig %>%
      filter(log2FoldChange <= -log2fc_cutoff)

    message("Upregulated:            ", nrow(up_genes))
    message("Downregulated:          ", nrow(down_genes))

    up_genes_csv <- file.path(RESULTS_DIR, out_name("ora_genes_up.csv"))
    down_genes_csv <- file.path(RESULTS_DIR, out_name("ora_genes_down.csv"))
    write_csv(select(up_genes, all_of(gene_cols)), up_genes_csv)
    write_csv(select(down_genes, all_of(gene_cols)), down_genes_csv)
    gene_csvs <- c(up_genes_csv, down_genes_csv)

    enrich_up <- run_enricher(query_ids(up_genes), background, sets, "up")
    enrich_down <- run_enricher(query_ids(down_genes), background, sets, "down")

    result_table <- bind_rows(
      enrich_to_table(enrich_up, "up", sets$ont_map, sets$id_type, gene_map),
      enrich_to_table(enrich_down, "down", sets$ont_map, sets$id_type, gene_map)
    )
    enrich_res <- list(up = enrich_up, down = enrich_down, params = params)
  } else {
    all_genes_csv <- file.path(RESULTS_DIR, out_name("ora_genes_all.csv"))
    write_csv(select(deg_sig, all_of(gene_cols)), all_genes_csv)
    gene_csvs <- all_genes_csv

    enrich_all <- run_enricher(query_ids(deg_sig), background, sets, "none")
    result_table <- enrich_to_table(enrich_all, "none", sets$ont_map, sets$id_type, gene_map)
    enrich_res <- list(none = enrich_all, params = params)
  }

  if (is.null(result_table) || nrow(result_table) == 0) {
    result_table <- empty_ora_table()
  }

  if (split_direction == "yes") {
    term_counts <- c(
      up = n_sig_terms(result_table %>% filter(.data$direction == "up")),
      down = n_sig_terms(result_table %>% filter(.data$direction == "down"))
    )
  } else {
    term_counts <- c(none = n_sig_terms(result_table))
  }

  sig_table <- significant_ora_table(result_table)
  if (nrow(sig_table) == 0) {
    if (split_direction == "yes") {
      message("No significant terms found for either direction.")
    } else {
      message("No significant terms found.")
    }
  }

  write_csv(sig_table, ora_csv)
  write_csv(result_table, ora_all_csv)
  saveRDS(enrich_res, ora_rds)

  message("================")
  message("Terms (enricher BH FDR < ", ORA_FDR, ")")
  message("  all tested: ", nrow(result_table))
  message("  significant:")
  for (label in names(term_counts)) {
    message("    ", label, ": ", term_counts[[label]])
  }
  message("CSV files written:")
  message("  ", ora_csv)
  message("  ", ora_all_csv)
  for (gene_csv in gene_csvs) {
    message("  ", gene_csv)
  }
  message("RDS file written:")
  message("  ", ora_rds)
}

if (sys.nframe() == 0L) {
  source(file.path(dirname(dirname(normalizePath(
    sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)),
    winslash = "/", mustWork = TRUE
  ))), "cli.R"))
  invoke_module(ora_enricher_options, ora_enricher_run)
}
