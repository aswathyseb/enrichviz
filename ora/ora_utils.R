# Shared helpers for ORA tables from gprofiler, clusterProfiler, fgsea, etc.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
})

DIRECTION_ALIASES <- c(
  "up", "down", "none", "all",
  "upregulated", "downregulated",
  "up-regulated", "down-regulated",
  "positive", "negative",
  "activated", "suppressed",
  "both", "significant"
)

first_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[[1]]
}

normalize_ont <- function(x) {
  x <- toupper(as.character(x))
  x <- str_replace(x, "^GO[_:]", "")
  recode(
    x,
    "BP" = "BP", "GO:BP" = "BP", "GOBP" = "BP",
    "BIOLOGICAL_PROCESS" = "BP", "BIOLOGICAL PROCESS" = "BP",
    "CC" = "CC", "GO:CC" = "CC", "GOCC" = "CC",
    "CELLULAR_COMPONENT" = "CC", "CELLULAR COMPONENT" = "CC",
    "MF" = "MF", "GO:MF" = "MF", "GOMF" = "MF",
    "MOLECULAR_FUNCTION" = "MF", "MOLECULAR FUNCTION" = "MF",
    .default = x
  )
}

normalize_direction <- function(x) {
  x <- tolower(as.character(x))
  recode(
    x,
    "up" = "up", "upregulated" = "up", "up-regulated" = "up",
    "positive" = "up", "activated" = "up",
    "down" = "down", "downregulated" = "down", "down-regulated" = "down",
    "negative" = "down", "suppressed" = "down",
    "none" = "none", "all" = "none", "both" = "none", "significant" = "none",
    .default = x
  )
}

looks_like_direction <- function(x) {
  vals <- unique(tolower(as.character(x)))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  length(vals) > 0 && all(vals %in% DIRECTION_ALIASES)
}

parse_ratio <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(p) {
    nums <- suppressWarnings(as.numeric(p))
    if (length(nums) >= 2 && all(is.finite(nums[1:2])) && nums[[2]] != 0) {
      nums[[1]] / nums[[2]]
    } else if (length(nums) >= 1) {
      nums[[1]]
    } else {
      NA_real_
    }
  }, numeric(1))
}

ratio_numerator <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(p) {
    nums <- suppressWarnings(as.numeric(p))
    if (length(nums) >= 1 && is.finite(nums[[1]])) nums[[1]] else NA_real_
  }, numeric(1))
}

ratio_denominator <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(p) {
    nums <- suppressWarnings(as.numeric(p))
    if (length(nums) >= 2 && is.finite(nums[[2]])) nums[[2]] else NA_real_
  }, numeric(1))
}

as_gene_string <- function(x) {
  if (is.list(x)) {
    return(vapply(x, function(v) paste(v, collapse = "/"), character(1)))
  }
  gsub(",", "/", as.character(x), fixed = TRUE)
}

has_slot <- function(x, slot) {
  isS4(x) && slot %in% methods::slotNames(x)
}

# Pull a result data.frame out of CSV-like tables, RDS objects, or S4 results.
extract_ora_df <- function(x) {
  if (is.data.frame(x)) {
    return(x)
  }

  if (has_slot(x, "compareClusterResult")) {
    return(as.data.frame(x@compareClusterResult))
  }
  if (has_slot(x, "result")) {
    return(as.data.frame(x@result))
  }

  if (is.list(x)) {
    wrap_names <- intersect(c("up", "down", "all", "none"), names(x))
    wrap_names <- wrap_names[vapply(x[wrap_names], function(v) !is.null(v), logical(1))]
    if (length(wrap_names) > 0 && is.null(x$result)) {
      pieces <- lapply(wrap_names, function(nm) {
        piece <- extract_ora_df(x[[nm]])
        if (!"direction" %in% colnames(piece)) {
          piece$direction <- nm
        }
        piece
      })
      return(bind_rows(pieces))
    }
    if (is.data.frame(x$result)) {
      return(x$result)
    }
    dfs <- x[vapply(x, is.data.frame, logical(1))]
    if (length(dfs) > 0) {
      return(bind_rows(dfs))
    }
  }

  stop("Could not extract an ORA result table from the input.", call. = FALSE)
}

read_ora_input <- function(path) {
  if (is.data.frame(path) || is.list(path) || isS4(path)) {
    return(extract_ora_df(path))
  }
  if (!is.character(path) || length(path) != 1) {
    stop("Input must be a file path or an ORA result object.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }

  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    return(extract_ora_df(readRDS(path)))
  }
  if (grepl("\\.(tsv|txt)$", path, ignore.case = TRUE)) {
    return(read_tsv(path, show_col_types = FALSE))
  }
  read_csv(path, show_col_types = FALSE)
}

infer_ora_format <- function(df) {
  nms <- colnames(df)
  if (all(c("term_id", "source") %in% nms) || all(c("term_id", "term_name") %in% nms)) {
    return("gprofiler")
  }
  if (all(c("ID", "Description") %in% nms) && any(c("p.adjust", "pvalue", "GeneRatio") %in% nms)) {
    return("clusterProfiler")
  }
  if (all(c("pathway", "padj") %in% nms) || all(c("pathway", "pval") %in% nms)) {
    return("fgsea")
  }
  "generic"
}

is_go_term_id <- function(x) {
  grepl("^GO:[0-9]+$", as.character(x), ignore.case = TRUE)
}

is_kegg_term_id <- function(x) {
  x <- sub("^path:", "", as.character(x), ignore.case = TRUE)
  grepl("^[a-z]{2,4}[0-9]{5}$", x, ignore.case = TRUE)
}

# Infer the gene-set database from term IDs / ontology columns (GO, KEGG, or generic ORA).
infer_ora_source <- function(df) {
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
  if (length(ids) > 0 && mean(is_kegg_term_id(ids)) >= 0.5) {
    return("KEGG")
  }
  "ORA"
}

go_ontologies <- function(ont) {
  ont <- unique(toupper(as.character(ont)))
  ont[ont %in% c("BP", "CC", "MF")]
}

# File-name stem and plot title from the detected database and optional GO ontology.
ora_plot_meta <- function(df, ont = NULL) {
  src <- infer_ora_source(df)
  ont <- go_ontologies(ont)
  if (src == "GO" && length(ont) == 1) {
    return(list(
      source = src,
      tag = paste0("ora_", ont),
      title = paste("GO", ont, "ORA"),
      term_label = paste("GO", ont)
    ))
  }
  if (src == "GO") {
    return(list(source = src, tag = "ora_GO", title = "GO ORA", term_label = "GO"))
  }
  if (src == "KEGG") {
    return(list(source = src, tag = "ora_KEGG", title = "KEGG ORA", term_label = "KEGG"))
  }
  list(source = src, tag = "ora", title = "ORA", term_label = "ORA")
}

ora_has_signed_direction <- function(plot_df) {
  all(c("up", "down") %in% unique(as.character(plot_df$direction)))
}

# When both up and down are present, leave room for negative values.
# Otherwise pin 0 to the start (or end, if all values are negative).
ora_x_scale <- function(plot_df, x) {
  if (ora_has_signed_direction(plot_df)) {
    return(scale_x_continuous())
  }
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0 || min(x) >= 0) {
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
  } else if (max(x) <= 0) {
    scale_x_continuous(limits = c(NA, 0), expand = expansion(mult = c(0.05, 0)))
  } else {
    scale_x_continuous()
  }
}

ora_zero_vline <- function(plot_df) {
  if (ora_has_signed_direction(plot_df)) {
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4)
  }
}

pick_direction_column <- function(df) {
  dir_col <- first_col(df, "direction")
  if (!is.na(dir_col)) {
    return(dir_col)
  }

  cand <- first_col(df, c("Cluster", "cluster"))
  if (!is.na(cand) && looks_like_direction(df[[cand]])) {
    return(cand)
  }

  NA_character_
}

# Map any supported ORA table onto a common set of columns.
normalize_ora_table <- function(df) {
  if (!is.data.frame(df)) {
    df <- extract_ora_df(df)
  }
  if (nrow(df) == 0) {
    return(tibble(
      ID = character(),
      Description = character(),
      p.adjust = double(),
      pvalue = double(),
      Count = double(),
      setSize = double(),
      querySize = double(),
      universeSize = double(),
      GeneRatio = double(),
      FoldEnrichment = double(),
      zScore = double(),
      ONTOLOGY = character(),
      direction = character(),
      geneID = character()
    ))
  }

  fmt <- infer_ora_format(df)
  message("Detected ORA table format: ", fmt)

  id_col <- first_col(df, c("term_id", "ID", "GO_ID", "go_id", "gs_exact_source", "pathway"))
  desc_col <- first_col(df, c("term_name", "Description", "description", "gs_name", "term"))
  adj_col <- first_col(df, c(
    "p.adjust", "p_adjust", "adjusted_p_value", "adjusted.p.value",
    "FDR", "padj", "qvalue", "q.value"
  ))
  raw_col <- first_col(df, c(
    "pvalue", "p_value", "pval", "p.value", "native_p_value", "native.p.value"
  ))
  ont_col <- first_col(df, c("ONTOLOGY", "ontology", "source", "ont"))
  dir_col <- pick_direction_column(df)
  count_col <- first_col(df, c("intersection_size", "Count", "count"))
  size_col <- first_col(df, c("term_size", "setSize", "set_size"))
  query_col <- first_col(df, c("query_size", "querySize"))
  universe_col <- first_col(df, c("effective_domain_size", "universeSize", "universe"))
  ratio_col <- first_col(df, c("GeneRatio", "gene_ratio"))
  bg_col <- first_col(df, c("BgRatio", "bg_ratio"))
  fe_col <- first_col(df, c("FoldEnrichment", "fold_enrichment", "foldEnrichment"))
  z_col <- first_col(df, c("zScore", "zscore", "z.score"))
  gene_col <- first_col(df, c(
    "intersection_genes", "geneID", "gene_id", "intersections", "intersection",
    "leadingEdge", "core_enrichment"
  ))

  if (is.na(id_col) && is.na(desc_col)) {
    stop(
      "Input must have a term ID or description column ",
      "(e.g. term_id / term_name, ID / Description, or pathway).",
      call. = FALSE
    )
  }
  if (is.na(adj_col) && is.na(raw_col)) {
    stop(
      "Input must have a p-value column ",
      "(e.g. p_value, adjusted_p_value, p.adjust, padj, or pvalue).",
      call. = FALSE
    )
  }

  ids <- if (!is.na(id_col)) as.character(df[[id_col]]) else as.character(df[[desc_col]])
  descs <- if (!is.na(desc_col)) as.character(df[[desc_col]]) else ids
  p_adj <- as.numeric(df[[if (!is.na(adj_col)) adj_col else raw_col]])
  p_raw <- as.numeric(df[[if (!is.na(raw_col)) raw_col else adj_col]])

  count <- if (!is.na(count_col)) as.numeric(df[[count_col]]) else NA_real_
  set_size <- if (!is.na(size_col)) as.numeric(df[[size_col]]) else NA_real_
  query_size <- if (!is.na(query_col)) as.numeric(df[[query_col]]) else NA_real_
  universe_size <- if (!is.na(universe_col)) as.numeric(df[[universe_col]]) else NA_real_
  gene_ratio <- if (!is.na(ratio_col)) parse_ratio(df[[ratio_col]]) else NA_real_
  fold_enrichment <- if (!is.na(fe_col)) as.numeric(df[[fe_col]]) else NA_real_
  z_score <- if (!is.na(z_col)) as.numeric(df[[z_col]]) else NA_real_

  if (all(is.na(count)) && !is.na(ratio_col)) {
    count <- ratio_numerator(df[[ratio_col]])
  }
  if (all(is.na(query_size)) && !is.na(ratio_col)) {
    query_size <- ratio_denominator(df[[ratio_col]])
  }
  if (all(is.na(set_size)) && !is.na(bg_col)) {
    set_size <- ratio_numerator(df[[bg_col]])
  }
  if (all(is.na(universe_size)) && !is.na(bg_col)) {
    universe_size <- ratio_denominator(df[[bg_col]])
  }
  if (all(is.na(count))) {
    count <- ifelse(is.na(set_size), 1, set_size)
  }
  if (all(is.na(set_size))) {
    set_size <- count
  }
  if (all(is.na(gene_ratio))) {
    gene_ratio <- ifelse(
      !is.na(query_size) & query_size > 0,
      count / query_size,
      NA_real_
    )
  }
  bg_ratio <- ifelse(
    !is.na(universe_size) & universe_size > 0,
    set_size / universe_size,
    NA_real_
  )
  if (all(is.na(fold_enrichment))) {
    fold_enrichment <- ifelse(
      is.finite(gene_ratio) & is.finite(bg_ratio) & bg_ratio > 0,
      gene_ratio / bg_ratio,
      NA_real_
    )
  }

  ont <- if (!is.na(ont_col)) normalize_ont(df[[ont_col]]) else NA_character_
  direction <- if (!is.na(dir_col)) normalize_direction(df[[dir_col]]) else "none"
  direction[is.na(direction) | direction == ""] <- "none"

  gene_id <- if (!is.na(gene_col)) as_gene_string(df[[gene_col]]) else NA_character_

  tibble(
    ID = ids,
    Description = ifelse(is.na(descs) | descs == "", ids, descs),
    p.adjust = p_adj,
    pvalue = p_raw,
    Count = count,
    setSize = set_size,
    querySize = query_size,
    universeSize = universe_size,
    GeneRatio = gene_ratio,
    FoldEnrichment = fold_enrichment,
    zScore = z_score,
    ONTOLOGY = ont,
    direction = direction,
    geneID = gene_id
  ) %>%
    filter(!is.na(ID), ID != "")
}

filter_ora_table <- function(df, ont = NULL, fdr_cutoff = 0.05) {
  ont <- go_ontologies(ont)
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
      )
    )
}

read_ora_table <- function(path) {
  normalize_ora_table(read_ora_input(path))
}

# Read any supported ORA result and keep significant terms.
# GO ontology filtering is applied only when BP/CC/MF is present in the table.
read_ora_plot_table <- function(path, ont = NULL, fdr_cutoff = 0.05) {
  df <- read_ora_table(path)
  src <- infer_ora_source(df)
  ont_filter <- if (src == "GO") go_ontologies(ont) else NULL
  if (src != "GO" && length(go_ontologies(ont)) == 1) {
    message("--ont ", go_ontologies(ont), " is GO-specific; plotting all ", src, " terms.")
  }
  filter_ora_table(df, ont = ont_filter, fdr_cutoff = fdr_cutoff)
}

ora_out_name <- function(stem, prefix = "") {
  prefix <- if (is.null(prefix)) "" else as.character(prefix)
  if (!nzchar(prefix)) {
    return(stem)
  }
  prefix <- gsub("[/\\\\]", "", prefix)
  sep <- if (grepl("[_-]$", prefix)) "" else "_"
  paste0(prefix, sep, stem)
}

# Load an ORA table, infer GO/KEGG, and apply FDR / GO-ontology filters.
prepare_ora_for_plot <- function(path, ont = NULL, fdr_cutoff = 0.05) {
  ora_all <- read_ora_table(path)
  plot_meta <- ora_plot_meta(ora_all, ont = ont)
  ont_filter <- if (plot_meta$source == "GO") go_ontologies(ont) else NULL
  if (plot_meta$source != "GO" && length(go_ontologies(ont)) == 1) {
    message("--ont ", go_ontologies(ont), " is GO-specific; plotting all ", plot_meta$source, " terms.")
  }
  ora_df <- filter_ora_table(ora_all, ont = ont_filter, fdr_cutoff = fdr_cutoff)
  list(df = ora_df, meta = plot_meta, ont_filter = ont_filter)
}
