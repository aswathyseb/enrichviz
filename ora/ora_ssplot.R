# Usage:
#   Rscript ora_ssplot.R --in gprofiler_GO.csv --outdir res -n 30 --ont BP --direction up
#   Rscript ora_ssplot.R --in kegg.csv --outdir res -n 30 --direction none [--prefix NAME]
#
# Semantic space plot of ORA terms (clusterProfiler/enrichplot ssplot-style).
# Uses Jaccard gene-overlap similarity (ssplot pairwise_termsim default) and
# projects terms into 2D with MDS. --direction selects which terms to plot:
# up, down, or none (top N overall).

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--in", type = "character", default = "gprofiler_GO.csv", dest = "input_file",
              metavar = "FILE",
              help = "ORA result CSV/TSV or RDS (gprofiler, clusterProfiler, fgsea)"),
  make_option("--outdir", type = "character", default = "res",
              metavar = "DIR", help = "output directory for PDF plots and tables"),
  make_option(c("-n", "--n"), type = "integer", default = 30, dest = "top_n",
              metavar = "N", help = "top terms to project [default: %default]"),
  make_option("--ont", type = "character", default = "BP",
              metavar = "BP|CC|MF|KEGG|all",
              help = paste(
                "GO ontology to plot, or KEGG/all to skip GO filtering",
                "[default: %default]"
              )),
  make_option("--direction", type = "character", default = "none",
              metavar = "up|down|none",
              help = paste(
                "term set to cluster: up or down = that direction after filtering,",
                "none = top N terms regardless of direction [default: %default]"
              )),
  make_option("--nCluster", type = "integer", default = NA, dest = "n_cluster",
              metavar = "N",
              help = "number of term clusters; default is floor(sqrt(n))"),
  make_option("--min-edge", type = "double", default = 0.2, dest = "min_edge",
              metavar = "NUM",
              help = "minimum Jaccard similarity to draw an edge [default: %default]"),
  make_option("--prefix", type = "character", default = "",
              metavar = "NAME",
              help = "optional prefix for output file names [default: none]")
)

opt_parser <- OptionParser(
  usage = "\nRscript ora_ssplot.R --in FILE --outdir DIR -n 30 [--ont BP] [--direction up] [--prefix NAME]\n",
  option_list = option_list,
  description = paste(
    "Semantic space plot of ORA terms from GO, KEGG, or another gene-set database using MDS.",
    "Term similarity is Jaccard overlap of genes in each term, matching ssplot() defaults.\n",
    "For GO tables, --ont BP, CC, or MF selects one ontology.",
    "For KEGG and other non-GO tables, ontology filtering is skipped. \n",
    "Use --direction up or down to plot that subset; none plots the top N terms overall.",
    "--prefix is prepended to output file names.",
    "",
    "Examples:",
    "  Rscript ora_ssplot.R --in gprofiler_GO.csv --outdir res -n 30 --ont BP --direction up",
    "  Rscript ora_ssplot.R --in kegg.csv --outdir res -n 30 --direction none --prefix kegg",
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
direction <- tolower(opt$direction)
n_cluster <- opt$n_cluster
min_edge <- opt$min_edge
prefix <- if (is.null(opt$prefix)) "" else as.character(opt$prefix)
fdr_cutoff <- 0.05

# enrichplot ssplot / set_enrichplot_color defaults
SSPLOT_COLORS <- c("#e06663", "#327eba")
N_WORDS <- 4
LABEL_WRAP <- 30
SSPLOT_STOP_WORDS <- c(
  "the", "and", "for", "with", "via", "by", "to", "a", "an", "in", "of", "on", "at"
)
GO_LEAD_WORDS <- c("positive", "negative")
GO_MOD_WORDS <- c("regulation")
GO_TAIL_WORDS <- c("process", "activity", "pathway", "function")

if (!ONT %in% c("BP", "CC", "MF", "KEGG", "ALL")) {
  stop("--ont must be BP, CC, MF, KEGG, or all", call. = FALSE)
}
if (!direction %in% c("up", "down", "none")) {
  stop("--direction must be up, down, or none", call. = FALSE)
}
if (!is.finite(top_n) || top_n < 1) {
  stop("-n must be a positive integer", call. = FALSE)
}
if (!is.na(n_cluster) && (!is.finite(n_cluster) || n_cluster < 1)) {
  stop("--nCluster must be a positive integer", call. = FALSE)
}
if (!is.finite(min_edge) || min_edge < 0 || min_edge > 1) {
  stop("--min-edge must be between 0 and 1", call. = FALSE)
}
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(ggforce)
  library(tidydr)
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
    arrange(p.adjust) %>%
    slice_head(n = n) %>%
    mutate(label = str_to_sentence(Description))

  if (anyDuplicated(plot_df$label) > 0) {
    plot_df <- plot_df %>%
      mutate(label = paste0(label, " (", ID, ")"))
  }

  plot_df
}

parse_gene_set <- function(x) {
  genes <- unlist(strsplit(as.character(x), "[,/;]"))
  genes <- str_trim(genes)
  unique(genes[!is.na(genes) & nzchar(genes)])
}

jaccard_similarity <- function(gene_sets) {
  n <- length(gene_sets)
  ids <- names(gene_sets)
  mat <- matrix(0, n, n, dimnames = list(ids, ids))
  for (i in seq_len(n)) {
    mat[i, i] <- 1
    if (i == n) {
      next
    }
    a <- gene_sets[[i]]
    for (j in (i + 1):n) {
      b <- gene_sets[[j]]
      u <- length(union(a, b))
      s <- if (u == 0) 0 else length(intersect(a, b)) / u
      mat[i, j] <- mat[j, i] <- s
    }
  }
  mat
}

# Match enrichplot::get_wordcloud() cleaning of GO descriptions.
clean_ssplot_terms <- function(terms) {
  x <- tolower(as.character(terms))
  x <- gsub(" in ", " ", x)
  x <- gsub(" [0-9]+ ", " ", x)
  x <- gsub("^[0-9]+ ", "", x)
  x <- gsub(" [0-9]+$", "", x)
  x <- gsub(" [a-z] ", " ", x)
  x <- gsub("^[a-z] ", "", x)
  x <- gsub(" [a-z]$", "", x)
  x <- gsub(" / ", " ", x)
  x <- gsub(" and ", " ", x)
  x <- gsub(" of ", " ", x)
  x <- gsub(",", " ", x)
  x <- gsub(" - ", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

# Match enrichplot::calculate_word_positions(): order words by where they
# typically appear in the original term names.
order_words_by_position <- function(terms, words) {
  term_words <- strsplit(clean_ssplot_terms(terms), "\\s+")
  ranks <- vapply(words, function(word) {
    hits <- unlist(lapply(term_words, function(tw) {
      idx <- which(tw == word)
      if (length(idx) == 0) {
        numeric()
      } else {
        idx[[1]]
      }
    }))
    if (length(hits) == 0) Inf else mean(hits)
  }, numeric(1))
  words[order(ranks, words)]
}

# After ssplot's bag-of-words, snap to GO-like English:
# "positive/negative regulation of <content> process/activity".
# Used only for GO tables; KEGG and other names keep word-frequency order.
naturalize_go_phrase <- function(words, terms) {
  words <- unique(tolower(words))
  words <- words[nzchar(words)]
  if (length(words) == 0) {
    return(character())
  }

  ordered <- order_words_by_position(terms, words)
  lead <- ordered[ordered %in% GO_LEAD_WORDS]
  mod <- ordered[ordered %in% GO_MOD_WORDS]
  tail <- ordered[ordered %in% GO_TAIL_WORDS]
  middle <- ordered[!ordered %in% c(GO_LEAD_WORDS, GO_MOD_WORDS, GO_TAIL_WORDS)]
  restore_go_prepositions(c(lead, mod, middle, tail))
}

restore_go_prepositions <- function(words) {
  if (length(words) == 0) {
    return(words)
  }
  out <- character()
  for (i in seq_along(words)) {
    out <- c(out, words[[i]])
    if (i == length(words)) {
      next
    }
    nxt <- words[[i + 1]]
    if (words[[i]] == "regulation" && nxt != "of") {
      out <- c(out, "of")
    }
    if (words[[i]] == "response" && !nxt %in% c("to", "of")) {
      out <- c(out, "to")
    }
  }
  out
}

# ssplot-style cluster tag: top frequent words, then rearranged to read naturally.
cluster_keywords <- function(descriptions, n_words = N_WORDS, go_style = TRUE) {
  descriptions <- as.character(descriptions)
  descriptions <- descriptions[!is.na(descriptions) & nzchar(descriptions)]
  if (length(descriptions) == 0) {
    return(NA_character_)
  }

  cleaned <- clean_ssplot_terms(descriptions)
  all_words <- unlist(strsplit(cleaned, "\\s+"))
  all_words <- all_words[nzchar(all_words)]
  if (length(all_words) == 0) {
    return(descriptions[[1]])
  }

  word_freq <- sort(table(all_words), decreasing = TRUE)
  keep <- names(word_freq)[!tolower(names(word_freq)) %in% SSPLOT_STOP_WORDS]
  if (length(keep) == 0) {
    keep <- names(word_freq)
  }
  top_words <- head(keep, n_words)
  phrase <- if (isTRUE(go_style)) {
    naturalize_go_phrase(top_words, descriptions)
  } else {
    order_words_by_position(descriptions, top_words)
  }
  paste(phrase, collapse = " ")
}

# KEGG (and other non-GO) names are already short titles. Mixing frequent
# words across pathways produces labels like "alzheimer ataxia disease".
# Use the pathway closest to the cluster centroid instead.
representative_cluster_label <- function(descriptions, x, y, p.adjust) {
  descriptions <- as.character(descriptions)
  ok <- !is.na(descriptions) & nzchar(descriptions)
  if (!any(ok)) {
    return(NA_character_)
  }
  descriptions <- descriptions[ok]
  x <- as.numeric(x)[ok]
  y <- as.numeric(y)[ok]
  p.adjust <- as.numeric(p.adjust)[ok]
  if (length(descriptions) == 1) {
    return(descriptions[[1]])
  }

  d <- (x - mean(x, na.rm = TRUE))^2 + (y - mean(y, na.rm = TRUE))^2
  d[!is.finite(d)] <- Inf
  p.adjust[!is.finite(p.adjust)] <- Inf
  descriptions[[order(d, p.adjust)[[1]]]]
}

wrap_label <- function(x, width = LABEL_WRAP) {
  vapply(as.character(x), function(s) {
    paste(strwrap(s, width = width), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

similarity_to_distance <- function(sim) {
  if (!isSymmetric(sim)) {
    sim <- (sim + t(sim)) / 2
  }
  sim[is.na(sim)] <- 0
  sim <- pmin(pmax(sim, 0), 1)
  diag(sim) <- 1
  offdiag <- row(sim) != col(sim)
  sim[offdiag & sim >= 1] <- 1 - .Machine$double.eps
  stats::as.dist(1 - sim)
}

mds_coords <- function(sim) {
  dist_mat <- similarity_to_distance(sim)
  n <- attr(dist_mat, "Size")
  if (n < 2) {
    stop("Need at least 2 terms for MDS.", call. = FALSE)
  }

  k <- min(2L, n - 1L)
  fit <- stats::cmdscale(dist_mat, k = k, eig = TRUE)
  pts <- as.matrix(fit$points)
  if (ncol(pts) == 1) {
    pts <- cbind(pts, 0)
  }
  colnames(pts) <- c("x", "y")
  rownames(pts) <- labels(dist_mat)

  eig <- as.numeric(fit$eig)
  eig <- eig[is.finite(eig) & eig > 0]
  list(coords = pts, eigenvalue = eig)
}

build_edges <- function(sim, min_edge) {
  ids <- rownames(sim)
  if (length(ids) < 2) {
    return(tibble(term1 = character(), term2 = character(), similarity = double()))
  }

  idx <- which(upper.tri(sim), arr.ind = TRUE)
  tibble(
    term1 = ids[idx[, 1]],
    term2 = ids[idx[, 2]],
    similarity = sim[idx]
  ) %>%
    filter(is.finite(similarity), similarity >= min_edge)
}

cluster_can_ellipse <- function(x, y) {
  xy <- cbind(as.numeric(x), as.numeric(y))
  xy <- xy[stats::complete.cases(xy), , drop = FALSE]
  xy <- unique(xy)
  if (nrow(xy) < 3) {
    return(FALSE)
  }
  covm <- stats::cov(xy)
  if (any(!is.finite(covm))) {
    return(FALSE)
  }
  ev <- eigen(covm, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, 0)
  if (max(ev) <= 0) {
    return(FALSE)
  }
  (min(ev) / max(ev)) > 1e-8
}

ora_ssplot <- function(nodes, edges, xlab, ylab, title, direction) {
  ellipse_df <- nodes %>%
    group_by(cluster) %>%
    filter(cluster_can_ellipse(x, y)) %>%
    ungroup()

  p <- ggplot(nodes, aes(x = x, y = y))

  if (nrow(ellipse_df) > 0) {
    p <- p +
      geom_mark_ellipse(
        data = ellipse_df,
        aes(group = cluster, fill = cluster_label),
        alpha = 0.3,
        color = NA,
        show.legend = TRUE
      )
  }

  if (nrow(edges) > 0) {
    p <- p +
      geom_segment(
        data = edges,
        aes(x = x, y = y, xend = xend, yend = yend),
        color = "grey",
        linewidth = 0.5,
        alpha = 0.8,
        inherit.aes = FALSE
      )
  }

  label_df <- nodes %>%
    group_by(cluster, cluster_label) %>%
    summarise(x = mean(x), y = mean(y), .groups = "drop") %>%
    mutate(label = wrap_label(cluster_label))

  p <- p +
    geom_point(aes(color = p.adjust, size = Count), alpha = 0.9) +
    geom_text_repel(
      data = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      bg.color = "white",
      bg.r = 0.1,
      size = 5,
      max.overlaps = Inf
    ) +
    scale_color_gradient(
      low = SSPLOT_COLORS[[1]],
      high = SSPLOT_COLORS[[2]],
      transform = "log10",
      name = "p.adjust",
      guide = guide_colorbar(reverse = TRUE, order = 2)
    ) +
    scale_size(range = c(3, 8), name = "Count", guide = guide_legend(order = 1))

  if (nrow(ellipse_df) > 0) {
    p <- p + scale_fill_discrete(name = "groups")
  }

  p +
    coord_equal() +
    labs(
      title = title,
      subtitle = if (direction == "none") {
        "MDS of Jaccard gene-overlap similarity"
      } else {
        paste("Direction:", direction, "| MDS of Jaccard gene-overlap similarity")
      },
      x = xlab,
      y = ylab
    ) +
    theme_dr()
}

message("Input:     ", input_file)
message("Output:    ", RESULTS_DIR)
message("Prefix:    ", if (nzchar(prefix)) prefix else "(none)")
message("N:         ", top_n)
message("Direction: ", direction)

loaded <- prepare_ora_for_plot(input_file, ont = ONT, fdr_cutoff = fdr_cutoff)
plot_meta <- loaded$meta
ont_filter <- loaded$ont_filter
ora_df <- loaded$df

message("Source:    ", plot_meta$source)
message("ONT:       ", if (length(ont_filter) == 1) ont_filter else plot_meta$source)

if (nrow(ora_df) == 0) {
  stop(
    "No ", plot_meta$term_label, " terms with FDR < ", fdr_cutoff, " in ", input_file,
    call. = FALSE
  )
}

plot_df <- select_top_terms(ora_df, top_n, direction)
if (nrow(plot_df) == 0) {
  if (direction %in% c("up", "down")) {
    present <- paste(unique(ora_df$direction), collapse = ", ")
    stop(
      "No ", direction, " ", plot_meta$term_label, " terms after filtering. ",
      "Directions present: ", present, ".",
      call. = FALSE
    )
  }
  stop("No terms left to plot after selecting the top ", top_n, ".", call. = FALSE)
}

plot_df <- plot_df %>%
  mutate(gene_list = lapply(geneID, parse_gene_set)) %>%
  filter(lengths(gene_list) > 0) %>%
  distinct(ID, .keep_all = TRUE)

if (nrow(plot_df) < 3) {
  stop(
    "Need at least 3 terms with gene lists for a Jaccard semantic space plot.",
    call. = FALSE
  )
}

message("Terms for MDS: ", nrow(plot_df))

gene_sets <- plot_df$gene_list
names(gene_sets) <- plot_df$ID
mat <- jaccard_similarity(gene_sets)

mds <- mds_coords(mat)
coords <- mds$coords[plot_df$ID, , drop = FALSE]

n_terms <- nrow(plot_df)
k <- if (is.na(n_cluster)) {
  max(1L, floor(sqrt(n_terms)))
} else {
  min(as.integer(n_cluster), n_terms)
}

set.seed(123)
cluster_id <- if (k == 1L) {
  rep(1L, n_terms)
} else {
  stats::kmeans(coords, centers = k)$cluster
}

nodes <- plot_df %>%
  mutate(
    x = coords[, "x"],
    y = coords[, "y"],
    cluster = as.integer(cluster_id)
  ) %>%
  group_by(cluster) %>%
  mutate(
    cluster_label = if (plot_meta$source == "GO") {
      cluster_keywords(Description, N_WORDS, go_style = TRUE)
    } else {
      representative_cluster_label(Description, x, y, p.adjust)
    }
  ) %>%
  ungroup()

cluster_names <- nodes %>% distinct(cluster, cluster_label)
if (anyDuplicated(cluster_names$cluster_label) > 0) {
  nodes <- nodes %>%
    mutate(cluster_label = paste0(cluster_label, " (", cluster, ")"))
}

edges <- build_edges(mat, min_edge)
if (nrow(edges) > 0) {
  xy <- nodes %>% transmute(ID, x, y)
  edges <- edges %>%
    left_join(xy, by = c("term1" = "ID")) %>%
    left_join(xy, by = c("term2" = "ID"), suffix = c("", "end")) %>%
    filter(is.finite(x), is.finite(y), is.finite(xend), is.finite(yend))
}

eig <- mds$eigenvalue
if (length(eig) >= 2 && sum(eig) > 0) {
  pct <- 100 * eig[seq_len(2)] / sum(eig)
  xlab <- paste0("Dimension1 (", format(pct[[1]], digits = 4), "%)")
  ylab <- paste0("Dimension2 (", format(pct[[2]], digits = 4), "%)")
} else {
  xlab <- "Dimension1"
  ylab <- "Dimension2"
}

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

if (direction == "none") {
  stem <- paste0(plot_meta$tag, "_ssplot")
} else {
  stem <- paste0(plot_meta$tag, "_", direction, "_ssplot")
}
plot_file <- file.path(RESULTS_DIR, ora_out_name(paste0(stem, ".pdf"), prefix))
plot_csv <- file.path(RESULTS_DIR, ora_out_name(paste0(stem, "_table.csv"), prefix))
sim_rds <- file.path(RESULTS_DIR, ora_out_name(paste0(stem, "_similarity.rds"), prefix))

write_csv(
  nodes %>%
    transmute(
      term_id = ID,
      term_name = Description,
      p.adjust,
      Count,
      direction,
      cluster,
      cluster_label,
      x,
      y
    ),
  plot_csv
)
saveRDS(mat, sim_rds)

p <- ora_ssplot(
  nodes, edges, xlab, ylab,
  title = paste(plot_meta$term_label, "semantic space"),
  direction = direction
)
ggsave(plot_file, plot = p, width = 12, height = 10)

dir_counts <- table(nodes$direction)

message("================")
message("Significant ", plot_meta$term_label, " terms (FDR < ", fdr_cutoff, "): ", nrow(ora_df))
message(
  "Plotted: ",
  paste(paste0(names(dir_counts), "=", as.integer(dir_counts)), collapse = ", ")
)
message("Clusters: ", n_distinct(nodes$cluster))
message("CSV files written:")
message("  ", plot_csv)
message("Plots created:")
message("  ", plot_file)
message("RDS file written:")
message("  ", sim_rds)
