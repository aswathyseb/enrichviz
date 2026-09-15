# Shared helpers for GSEA GO/KEGG scripts and plot tables.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(readr)
})

# 1. Create a ranked gene list 
# Named, decreasingly sorted gene list for clusterProfiler GSEA.
# Names are Ensembl IDs; values are the ranking statistic.
create_gene_list <- function(de, ranking_method = "signed_logp") {
  required_cols <- c("name", "gene", "log2FoldChange", "PValue")
  missing_cols <- setdiff(required_cols, colnames(de))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  de <- de %>%
    mutate(
      ensembl_gene = sub("\\..*$", "", name),
      gene_symbol = gene
    )

  if (ranking_method == "signed_logp") {
    # Protect against exact zero P-values
    min_nonzero_p <- min(de$PValue[de$PValue > 0], na.rm = TRUE)
    replacement_p <- min_nonzero_p / 10

    de <- de %>%
      mutate(
        PValue_for_rank = ifelse(PValue == 0, replacement_p, PValue),
        rank_score = sign(log2FoldChange) * -log10(PValue_for_rank)
      )
  } else if (ranking_method == "log2FC") {
    de <- de %>%
      mutate(rank_score = log2FoldChange)
  } else {
    stop("ranking_method must be 'signed_logp' or 'log2FC'")
  }

  rank_df <- de %>%
    filter(
      !is.na(ensembl_gene),
      ensembl_gene != "",
      !is.na(rank_score),
      is.finite(rank_score)
    ) %>%
    select(
      ensembl_gene,
      gene_symbol,
      log2FoldChange,
      PValue,
      FDR,
      rank_score
    ) %>%
    arrange(desc(abs(rank_score))) %>%
    distinct(ensembl_gene, .keep_all = TRUE) %>%
    arrange(desc(rank_score))

  gene_list <- rank_df$rank_score
  names(gene_list) <- rank_df$ensembl_gene
  gene_list <- sort(gene_list, decreasing = TRUE)

  if (anyDuplicated(names(gene_list)) > 0) {
    stop("Duplicated Ensembl IDs remain in the ranked gene list.")
  }

  list(gene_list = gene_list, rank_df = rank_df)
}

# 2. Parse the expression matrix from the differential expression result table
# Normalized expression matrix for GseaVis:
# first column = gene symbol, remaining columns = samples after falsePos.
parse_expr <- function(de, gene_col = "gene", after_col = "falsePos") {
  if (!gene_col %in% colnames(de)) {
    stop("Gene symbol column not found: ", gene_col)
  }
  if (!after_col %in% colnames(de)) {
    stop("Column not found: ", after_col)
  }

  after_idx <- match(after_col, colnames(de))
  if (after_idx >= ncol(de)) {
    stop("No expression columns found after '", after_col, "'")
  }

  expr_cols <- colnames(de)[(after_idx + 1):ncol(de)]

  expr <- de %>%
    transmute(
      gene_symbol = as.character(.data[[gene_col]]),
      across(all_of(expr_cols), as.numeric)
    ) %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    distinct(gene_symbol, .keep_all = TRUE)

  # GseaVis uses exp[,1]; that must be a character vector, not a tibble column.
  as.data.frame(expr, stringsAsFactors = FALSE)
}

# 3. Get highest positive-NES term ID from a clusterProfiler result table.
top_positive_nes_id <- function(result_table) {
  pos_nes <- result_table[result_table$NES > 0, , drop = FALSE]
  if (nrow(pos_nes) == 0) {
    stop("No terms with positive NES to plot.")
  }
  pos_nes$ID[which.max(pos_nes$NES)]
}

# 4. NES lollipop of top positive and negative terms (clusterProfiler result table).
nes_lollipop_plot <- function(result_table, ont, top_n = 10, fdr_cutoff = 0.05) {
  plot_df <- result_table %>%
    filter(!is.na(p.adjust), p.adjust < fdr_cutoff) %>%
    mutate(
      direction = ifelse(NES > 0, "up", "down"),
      label = str_to_sentence(Description)
    ) %>%
    group_by(direction) %>%
    slice_max(order_by = abs(NES), n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(NES) %>%
    mutate(label = factor(label, levels = label))

  if (nrow(plot_df) == 0) {
    stop("No significant terms to plot in the NES lollipop.")
  }

  ggplot(plot_df, aes(x = NES, y = label, color = p.adjust)) +
    geom_segment(aes(x = 0, xend = NES, yend = label), linewidth = 0.6, color = "grey70") +
    geom_point(aes(size = setSize), alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    scale_color_viridis_c(option = "magma", direction = -1, trans = "log10") +
    labs(
      title = paste("GO", ont, "GSEA"),
      subtitle = "Top positively and negatively enriched terms",
      x = "Normalized enrichment score (NES)",
      y = NULL,
      color = "p.adjust",
      size = "Gene-set size"
    ) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 10))
}

# ---------------------------------------------------------------------------
# GSEA result tables from fgsea or clusterProfiler (CSV/TSV/RDS).
# ---------------------------------------------------------------------------

gsea_first_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[[1]]
}

gsea_as_gene_string <- function(x) {
  if (is.list(x)) {
    return(vapply(x, function(v) paste(v, collapse = "/"), character(1)))
  }
  gsub(",", "/", as.character(x), fixed = TRUE)
}

gsea_leading_edge_count <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(p) {
    p <- trimws(p)
    p <- p[nzchar(p) & !is.na(p)]
    as.integer(length(p))
  }, integer(1))
}

gsea_has_slot <- function(x, slot) {
  isS4(x) && slot %in% methods::slotNames(x)
}

# Pull a result data.frame out of CSV-like tables, RDS objects, or S4 results.
extract_gsea_df <- function(x) {
  if (is.data.frame(x)) {
    return(x)
  }

  if (gsea_has_slot(x, "result")) {
    return(as.data.frame(x@result))
  }

  if (is.list(x)) {
    if (is.data.frame(x$table)) {
      return(x$table)
    }
    if (is.list(x$by_ontology)) {
      pieces <- lapply(x$by_ontology, function(v) {
        if (is.list(v) && is.data.frame(v$table)) {
          return(v$table)
        }
        extract_gsea_df(v)
      })
      return(bind_rows(pieces))
    }
    dfs <- x[vapply(x, is.data.frame, logical(1))]
    if (length(dfs) > 0) {
      return(bind_rows(dfs))
    }
  }

  stop("Could not extract a GSEA result table from the input.", call. = FALSE)
}

read_gsea_input <- function(path) {
  if (is.data.frame(path) || is.list(path) || isS4(path)) {
    return(extract_gsea_df(path))
  }
  if (!is.character(path) || length(path) != 1) {
    stop("Input must be a file path or a GSEA result object.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }

  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    return(extract_gsea_df(readRDS(path)))
  }
  if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) {
    return(read_tsv(path, show_col_types = FALSE))
  }
  read_csv(path, show_col_types = FALSE)
}

# Map fgsea or clusterProfiler GSEA tables onto a common set of columns.
normalize_gsea_table <- function(df) {
  if (!is.data.frame(df)) {
    df <- extract_gsea_df(df)
  }
  if (nrow(df) == 0) {
    return(tibble::tibble(
      ID = character(),
      Description = character(),
      p.adjust = double(),
      pvalue = double(),
      NES = double(),
      ES = double(),
      setSize = double(),
      ONTOLOGY = character(),
      geneID = character(),
      leadingEdgeCount = integer(),
      leadingEdgeRatio = double()
    ))
  }

  id_col <- gsea_first_col(df, c("ID", "pathway", "term_id", "GO_ID", "go_id"))
  desc_col <- gsea_first_col(df, c("Description", "term_name", "pathway", "gs_name"))
  adj_col <- gsea_first_col(df, c("p.adjust", "padj", "FDR", "qvalue"))
  raw_col <- gsea_first_col(df, c("pvalue", "pval", "p_value", "p.value"))
  nes_col <- gsea_first_col(df, "NES")
  es_col <- gsea_first_col(df, c("ES", "enrichmentScore"))
  size_col <- gsea_first_col(df, c("setSize", "size"))
  ont_col <- gsea_first_col(df, c("ONTOLOGY", "ontology"))
  gene_col <- gsea_first_col(df, c("leadingEdge", "core_enrichment", "geneID"))

  if (is.na(nes_col)) {
    stop("Input must have an NES column.", call. = FALSE)
  }
  if (is.na(adj_col) && is.na(raw_col)) {
    stop(
      "Input must have an adjusted or raw p-value column ",
      "(e.g. padj, p.adjust, FDR, pval, or pvalue).",
      call. = FALSE
    )
  }

  ids <- if (!is.na(id_col)) as.character(df[[id_col]]) else as.character(df[[desc_col]])
  descs <- if (!is.na(desc_col)) as.character(df[[desc_col]]) else ids
  descs <- ifelse(
    is.na(descs) | descs == "",
    ids,
    ifelse(grepl("_", descs, fixed = TRUE), str_to_sentence(gsub("_", " ", descs)), descs)
  )

  tibble::tibble(
    ID = ids,
    Description = descs,
    p.adjust = as.numeric(df[[if (!is.na(adj_col)) adj_col else raw_col]]),
    pvalue = as.numeric(df[[if (!is.na(raw_col)) raw_col else adj_col]]),
    NES = as.numeric(df[[nes_col]]),
    ES = if (!is.na(es_col)) as.numeric(df[[es_col]]) else NA_real_,
    setSize = if (!is.na(size_col)) as.numeric(df[[size_col]]) else NA_real_,
    ONTOLOGY = if (!is.na(ont_col)) toupper(as.character(df[[ont_col]])) else NA_character_,
    geneID = if (!is.na(gene_col)) gsea_as_gene_string(df[[gene_col]]) else NA_character_
  ) %>%
    mutate(
      leadingEdgeCount = ifelse(
        is.na(geneID) | geneID == "",
        NA_integer_,
        gsea_leading_edge_count(geneID)
      ),
      leadingEdgeRatio = ifelse(
        is.finite(setSize) & setSize > 0 & is.finite(leadingEdgeCount),
        leadingEdgeCount / setSize,
        NA_real_
      )
    ) %>%
    filter(!is.na(ID), ID != "", !is.na(NES), is.finite(NES))
}

is_go_term_id <- function(x) {
  grepl("^GO:[0-9]+$", as.character(x), ignore.case = TRUE)
}

is_kegg_term_id <- function(x) {
  x <- sub("^path:", "", as.character(x), ignore.case = TRUE)
  grepl("^[a-z]{2,4}[0-9]{5}$", x, ignore.case = TRUE) |
    grepl("^N[0-9]{5}$", x, ignore.case = TRUE)
}

# Infer GO vs KEGG from term IDs or the ONTOLOGY column.
infer_gsea_source <- function(df) {
  ids <- unique(as.character(df$ID))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  ont <- unique(na.omit(as.character(df$ONTOLOGY)))
  ont <- ont[nzchar(ont)]

  if (length(ids) > 0 && mean(is_go_term_id(ids)) >= 0.5) {
    return("GO")
  }
  if (length(ont) > 0 && all(toupper(ont) %in% c("BP", "CC", "MF"))) {
    return("GO")
  }
  if (length(ont) > 0 && all(toupper(ont) == "KEGG")) {
    return("KEGG")
  }
  if (length(ids) > 0 && mean(is_kegg_term_id(ids)) >= 0.5) {
    return("KEGG")
  }
  "GSEA"
}

gsea_go_ontologies <- function(ont) {
  ont <- unique(toupper(as.character(ont)))
  ont[ont %in% c("BP", "CC", "MF")]
}

# File-name stem and plot title from the detected database and optional GO ontology.
gsea_plot_meta <- function(df, ont = NULL) {
  src <- infer_gsea_source(df)
  ont <- gsea_go_ontologies(ont)
  if (src == "GO" && length(ont) == 1) {
    return(list(
      source = src,
      tag = paste0("gsea_", ont),
      title = paste("GO", ont, "GSEA"),
      term_label = paste("GO", ont)
    ))
  }
  if (src == "GO") {
    return(list(source = src, tag = "gsea_GO", title = "GO GSEA", term_label = "GO"))
  }
  if (src == "KEGG") {
    return(list(source = src, tag = "gsea_KEGG", title = "KEGG GSEA", term_label = "KEGG"))
  }
  list(source = src, tag = "gsea", title = "GSEA", term_label = "GSEA")
}

filter_gsea_table <- function(df, ont = NULL, fdr_cutoff = 0.05) {
  ont <- gsea_go_ontologies(ont)
  if (length(ont) == 1 && "ONTOLOGY" %in% colnames(df)) {
    present <- unique(na.omit(df$ONTOLOGY))
    if (ont %in% present) {
      df <- df %>% filter(ONTOLOGY == ont)
    } else if (length(present) > 0) {
      message(
        "Requested ontology ", ont,
        " was not found; keeping ontologies present in the table: ",
        paste(present, collapse = ", ")
      )
    }
  }

  df %>%
    filter(
      !is.na(p.adjust),
      is.finite(p.adjust),
      p.adjust < fdr_cutoff
    ) %>%
    mutate(
      p.adjust = ifelse(
        p.adjust <= 0,
        min(c(p.adjust[p.adjust > 0], .Machine$double.xmin), na.rm = TRUE) / 10,
        p.adjust
      ),
      direction = ifelse(NES > 0, "up", "down")
    )
}

read_gsea_table <- function(path) {
  df <- normalize_gsea_table(read_gsea_input(path))
  src <- infer_gsea_source(df)
  if (src == "KEGG") {
    df$ONTOLOGY <- ifelse(
      is.na(df$ONTOLOGY) | df$ONTOLOGY == "",
      "KEGG",
      df$ONTOLOGY
    )
  }
  df
}

# Read any supported GSEA result and keep significant terms.
# GO ontology filtering is applied only when BP/CC/MF is present in the table.
read_gsea_plot_table <- function(path, ont = NULL, fdr_cutoff = 0.05) {
  prepare_gsea_for_plot(path, ont = ont, fdr_cutoff = fdr_cutoff)$df
}

# Load a GSEA table, infer GO/KEGG, and apply FDR / GO-ontology filters.
prepare_gsea_for_plot <- function(path, ont = NULL, fdr_cutoff = 0.05) {
  gsea_all <- read_gsea_table(path)
  src <- infer_gsea_source(gsea_all)
  if (src == "GO" && all(is.na(gsea_all$ONTOLOGY) | gsea_all$ONTOLOGY == "")) {
    ont_stamp <- gsea_go_ontologies(ont)
    if (length(ont_stamp) == 1) {
      gsea_all$ONTOLOGY <- ont_stamp
    }
  }
  plot_meta <- gsea_plot_meta(gsea_all, ont = ont)
  ont_filter <- if (plot_meta$source == "GO") gsea_go_ontologies(ont) else NULL
  if (plot_meta$source != "GO" && length(gsea_go_ontologies(ont)) == 1) {
    message(
      "--ont ", gsea_go_ontologies(ont),
      " is GO-specific; plotting all ", plot_meta$source, " terms."
    )
  }
  gsea_df <- filter_gsea_table(gsea_all, ont = ont_filter, fdr_cutoff = fdr_cutoff)
  list(df = gsea_df, all = gsea_all, meta = plot_meta, ont_filter = ont_filter)
}

# Named ranking statistic for ridgeplots: one row per gene ID (symbol or Ensembl).
read_gsea_rank_table <- function(path) {
  if (is.numeric(path) && !is.null(names(path))) {
    return(tibble::tibble(gene = names(path), rank_score = as.numeric(unname(path))))
  }
  if (is.data.frame(path)) {
    df <- path
  } else if (is.character(path) && length(path) == 1) {
    if (!file.exists(path)) {
      stop("Ranked-gene file not found: ", path, call. = FALSE)
    }
    if (grepl("\\.rds$", path, ignore.case = TRUE)) {
      obj <- readRDS(path)
      if (is.numeric(obj) && !is.null(names(obj))) {
        return(tibble::tibble(gene = names(obj), rank_score = as.numeric(unname(obj))))
      }
      if (is.list(obj) && is.numeric(obj$stats) && !is.null(names(obj$stats))) {
        return(tibble::tibble(gene = names(obj$stats), rank_score = as.numeric(unname(obj$stats))))
      }
      if (is.list(obj) && is.numeric(obj$gene_list) && !is.null(names(obj$gene_list))) {
        return(tibble::tibble(
          gene = names(obj$gene_list),
          rank_score = as.numeric(unname(obj$gene_list))
        ))
      }
      if (is.data.frame(obj)) {
        df <- obj
      } else {
        stop("Could not extract a ranked gene list from: ", path, call. = FALSE)
      }
    } else if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) {
      df <- read_tsv(path, show_col_types = FALSE)
    } else {
      df <- read_csv(path, show_col_types = FALSE)
    }
  } else {
    stop("Ranked genes must be a file path or a named numeric vector.", call. = FALSE)
  }

  score_col <- gsea_first_col(df, c("rank_score", "stat", "ranking_metric", "log2FoldChange"))
  if (is.na(score_col)) {
    stop(
      "Ranked-gene table must have a score column ",
      "(e.g. rank_score, stat, or log2FoldChange).",
      call. = FALSE
    )
  }

  scores <- as.numeric(df[[score_col]])
  pieces <- list()
  sym_col <- gsea_first_col(df, c("gene_symbol", "gene", "symbol", "SYMBOL"))
  ens_col <- gsea_first_col(df, c("ensembl_gene", "ensembl", "ENSEMBL", "name"))
  if (!is.na(sym_col)) {
    pieces[[length(pieces) + 1]] <- tibble::tibble(
      gene = as.character(df[[sym_col]]),
      rank_score = scores
    )
  }
  if (!is.na(ens_col)) {
    pieces[[length(pieces) + 1]] <- tibble::tibble(
      gene = sub("\\..*$", "", as.character(df[[ens_col]])),
      rank_score = scores
    )
  }
  if (length(pieces) == 0) {
    stop(
      "Ranked-gene table must have a gene column ",
      "(e.g. gene_symbol or ensembl_gene).",
      call. = FALSE
    )
  }

  bind_rows(pieces) %>%
    filter(!is.na(gene), gene != "", is.finite(rank_score)) %>%
    distinct(gene, .keep_all = TRUE)
}

# Expand leading-edge genes and join them to ranking scores.
expand_gsea_leading_edge <- function(df, ranks) {
  if (!"geneID" %in% colnames(df) || all(is.na(df$geneID) | df$geneID == "")) {
    stop("Leading-edge genes (geneID) are required for a ridgeplot.", call. = FALSE)
  }
  if (!all(c("gene", "rank_score") %in% colnames(ranks))) {
    stop("Rank table must have gene and rank_score columns.", call. = FALSE)
  }

  pieces <- lapply(seq_len(nrow(df)), function(i) {
    genes <- trimws(unlist(strsplit(as.character(df$geneID[[i]]), "/", fixed = TRUE)))
    genes <- genes[nzchar(genes) & !is.na(genes)]
    if (length(genes) == 0) {
      return(NULL)
    }
    out <- df[rep(i, length(genes)), , drop = FALSE]
    out$gene <- genes
    out
  })
  gene_df <- bind_rows(pieces)
  if (nrow(gene_df) == 0) {
    stop("No leading-edge genes were found in the selected terms.", call. = FALSE)
  }

  ranks <- ranks %>%
    mutate(gene_key = toupper(gene)) %>%
    filter(!is.na(gene_key), gene_key != "", is.finite(rank_score)) %>%
    distinct(gene_key, .keep_all = TRUE) %>%
    select(gene_key, rank_score)

  gene_df %>%
    mutate(gene_key = toupper(gene)) %>%
    inner_join(ranks, by = "gene_key")
}

normalize_go_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("_", ":", x, fixed = TRUE)
  x <- sub("^GO:+", "GO:", x)
  ifelse(grepl("^GO:", x), x, paste0("GO:", x))
}

# Keep GO IDs canonical; leave KEGG pathway IDs (hsa03010, N01394) unchanged.
normalize_gsea_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^path:", "", x, ignore.case = TRUE)
  x <- gsub("_", ":", x, fixed = TRUE)
  goish <- grepl("^GO:", x, ignore.case = TRUE) | grepl("^[0-9]{7}$", x)
  ifelse(goish, normalize_go_id(x), x)
}

gsea_ids_match <- function(ids, query) {
  toupper(normalize_gsea_id(ids)) == toupper(normalize_gsea_id(query))
}

gsea_term_tag <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", normalize_gsea_id(x))
}

split_gene_ids <- function(x) {
  genes <- trimws(unlist(strsplit(as.character(x), "/", fixed = TRUE)))
  genes[nzchar(genes) & !is.na(genes)]
}

# Weighted running enrichment score (fgsea / clusterProfiler, exponent 1).
gsea_running_es <- function(stats, pathway_genes) {
  if (is.null(stats) || length(stats) == 0 || is.null(names(stats))) {
    stop("A named ranked gene list is required to draw the running ES curve.", call. = FALSE)
  }
  stats <- sort(stats, decreasing = TRUE)
  pathway_genes <- unique(as.character(pathway_genes))
  pathway_genes <- pathway_genes[nzchar(pathway_genes) & !is.na(pathway_genes)]

  in_set <- names(stats) %in% pathway_genes
  n_hit <- sum(in_set)
  if (n_hit == 0) {
    stop("No genes from the selected term overlap the ranked list.", call. = FALSE)
  }

  n <- length(stats)
  n_miss <- n - n_hit
  abs_stat <- abs(as.numeric(stats))
  nr <- sum(abs_stat[in_set])
  if (!is.finite(nr) || nr <= 0) {
    nr <- n_hit
  }
  step <- ifelse(in_set, abs_stat / nr, if (n_miss > 0) -1 / n_miss else 0)

  tibble::tibble(
    rank = seq_len(n),
    gene = names(stats),
    stat = as.numeric(stats),
    in_set = in_set,
    running_es = cumsum(step)
  )
}

gsea_named_stats <- function(x) {
  if (is.numeric(x) && !is.null(names(x))) {
    return(sort(x, decreasing = TRUE))
  }
  NULL
}

gsea_align_ids <- function(ids, stats, rank_file = NULL) {
  ids <- unique(as.character(ids))
  ids <- ids[nzchar(ids) & !is.na(ids)]
  if (length(ids) == 0 || is.null(stats)) {
    return(ids)
  }
  if (sum(ids %in% names(stats)) > 0) {
    return(ids)
  }

  nm <- names(stats)
  hit <- match(toupper(ids), toupper(nm))
  mapped <- nm[hit]
  mapped <- mapped[!is.na(mapped)]
  if (length(mapped) > 0) {
    return(unique(mapped))
  }

  if (is.null(rank_file) || !file.exists(rank_file) ||
      grepl("\\.rds$", rank_file, ignore.case = TRUE)) {
    return(ids)
  }

  rank_df <- read_csv(rank_file, show_col_types = FALSE)
  sym_col <- gsea_first_col(rank_df, c("gene_symbol", "gene", "symbol", "SYMBOL"))
  ens_col <- gsea_first_col(rank_df, c("ensembl_gene", "ensembl", "ENSEMBL", "name"))
  if (is.na(sym_col) || is.na(ens_col)) {
    return(ids)
  }

  ens <- sub("\\..*$", "", as.character(rank_df[[ens_col]]))
  sym <- as.character(rank_df[[sym_col]])
  to_ens <- setNames(ens, toupper(sym))
  mapped <- unname(to_ens[toupper(ids)])
  mapped <- mapped[!is.na(mapped) & mapped != ""]
  if (sum(mapped %in% names(stats)) > 0) {
    return(unique(mapped))
  }
  to_sym <- setNames(sym, toupper(ens))
  mapped <- unname(to_sym[toupper(ids)])
  mapped <- mapped[!is.na(mapped) & mapped != ""]
  if (sum(mapped %in% names(stats)) > 0) {
    return(unique(mapped))
  }
  ids
}

# Stats, full gene set, and leading-edge genes for one GO or KEGG term.
load_gsea_nes_term <- function(path, go_id, rank_file = NULL) {
  go_id <- normalize_gsea_id(go_id)
  rds_path <- NULL
  obj <- NULL

  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    rds_path <- path
  } else {
    guess <- sub("\\.(csv|tsv|txt)$", ".rds", path, ignore.case = TRUE)
    if (file.exists(guess)) {
      rds_path <- guess
    }
  }
  if (!is.null(rds_path) && file.exists(rds_path)) {
    obj <- readRDS(rds_path)
  }

  table <- if (is.null(obj)) {
    read_gsea_table(path)
  } else {
    df <- normalize_gsea_table(extract_gsea_df(obj))
    src <- infer_gsea_source(df)
    if (src == "KEGG") {
      df$ONTOLOGY <- ifelse(
        is.na(df$ONTOLOGY) | df$ONTOLOGY == "",
        "KEGG",
        df$ONTOLOGY
      )
    }
    df
  }
  hit <- which(gsea_ids_match(table$ID, go_id))
  if (length(hit) == 0) {
    stop("Term ID not found in GSEA results: ", go_id, call. = FALSE)
  }
  term <- table[hit[[1]], , drop = FALSE]

  stats <- NULL
  pathway_genes <- NULL

  if (gsea_has_slot(obj, "geneList")) {
    stats <- gsea_named_stats(obj@geneList)
  }
  if (gsea_has_slot(obj, "geneSets")) {
    gs_names <- names(obj@geneSets)
    gs_hit <- which(gsea_ids_match(gs_names, go_id))
    if (length(gs_hit) > 0) {
      pathway_genes <- obj@geneSets[[gs_hit[[1]]]]
    }
  }

  if (is.list(obj) && !isS4(obj)) {
    if (is.null(stats)) {
      stats <- gsea_named_stats(obj$stats)
    }
    pieces <- list()
    if (is.list(obj$by_ontology)) {
      pieces <- obj$by_ontology
    } else {
      pieces <- list(obj)
    }
    for (piece in pieces) {
      pws <- piece$pathways
      raw_tab <- if (is.data.frame(piece$table)) piece$table else NULL
      if (is.null(pws) || is.null(raw_tab)) {
        next
      }
      raw_ids <- if ("ID" %in% colnames(raw_tab)) raw_tab$ID else raw_tab[[1]]
      raw_hit <- which(gsea_ids_match(raw_ids, go_id))
      if (length(raw_hit) == 0) {
        next
      }
      row <- raw_tab[raw_hit[[1]], , drop = FALSE]
      keys <- unique(c(
        as.character(row$pathway),
        as.character(row$ID),
        as.character(row$gs_name)
      ))
      keys <- keys[keys %in% names(pws)]
      if (length(keys) > 0) {
        pathway_genes <- pws[[keys[[1]]]]
        break
      }
    }
    if (is.null(stats) && is.list(obj$pathways) && is.numeric(obj$stats)) {
      stats <- gsea_named_stats(obj$stats)
    }
  }

  if (is.null(stats) && !is.null(rank_file) && file.exists(rank_file)) {
    if (grepl("\\.rds$", rank_file, ignore.case = TRUE)) {
      rank_obj <- readRDS(rank_file)
      stats <- gsea_named_stats(rank_obj)
      if (is.null(stats) && is.list(rank_obj)) {
        stats <- gsea_named_stats(rank_obj$stats)
        if (is.null(stats)) {
          stats <- gsea_named_stats(rank_obj$gene_list)
        }
      }
    }
    if (is.null(stats)) {
      rank_df <- read_csv(rank_file, show_col_types = FALSE)
      ens_col <- gsea_first_col(rank_df, c("ensembl_gene", "ensembl", "name"))
      score_col <- gsea_first_col(rank_df, c("rank_score", "stat", "log2FoldChange"))
      if (!is.na(ens_col) && !is.na(score_col)) {
        stats <- setNames(
          as.numeric(rank_df[[score_col]]),
          sub("\\..*$", "", as.character(rank_df[[ens_col]]))
        )
        stats <- stats[is.finite(stats) & !is.na(names(stats)) & names(stats) != ""]
        stats <- tapply(stats, names(stats), `[`, 1)
        stats <- gsea_named_stats(stats)
      }
    }
  }

  leading_edge <- split_gene_ids(term$geneID)
  if (is.null(pathway_genes) || length(pathway_genes) == 0) {
    pathway_genes <- leading_edge
    message(
      "Full gene set was not found in the GSEA RDS; ",
      "the running ES curve uses leading-edge genes only. ",
      "Pass an fgsea or clusterProfiler RDS for the complete curve."
    )
  }
  pathway_genes <- gsea_align_ids(pathway_genes, stats, rank_file)
  leading_edge <- gsea_align_ids(leading_edge, stats, rank_file)

  list(
    term = term,
    go_id = go_id,
    stats = stats,
    pathway_genes = unique(as.character(pathway_genes)),
    leading_edge = leading_edge
  )
}

# Prepend an optional --prefix to an output file name.
gsea_out_name <- function(stem, prefix = "") {
  prefix <- if (is.null(prefix)) "" else as.character(prefix)
  if (!nzchar(prefix)) {
    return(stem)
  }
  prefix <- gsub("[/\\\\]", "", prefix)
  sep <- if (grepl("[_-]$", prefix)) "" else "_"
  paste0(prefix, sep, stem)
}
