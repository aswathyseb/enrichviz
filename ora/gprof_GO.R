# Usage:
#   Rscript gprof_GO.R --in edger.csv --outdir res [--fdr 0.05] [--log2fc 1] [--direction no] [--prefix NAME]
#
# Select significant DEGs by FDR and |log2FC|, then run gprofiler2 ORA
# across GO BP, MF, and CC. With --direction yes, up- and down-regulated
# genes are tested separately; with --direction no (default), all
# significant genes are tested together.
# --prefix is prepended to output file names; default is none.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--in", type = "character", default = "edger.csv", dest = "input_file",
              metavar = "FILE", help = "edgeR differential expression CSV"),
  make_option("--outdir", type = "character", default = "gprof_go",
              metavar = "DIR", help = "output directory for tables and RDS objects"),
  make_option("--fdr", type = "double", default = 0.05, dest = "fdr_cutoff",
              metavar = "NUM", help = "DEG FDR cutoff [default: %default]"),
  make_option("--log2fc", type = "double", default = 1, dest = "log2fc_cutoff",
              metavar = "NUM", help = "DEG |log2FoldChange| cutoff [default: %default]"),
  make_option("--organism", type = "character", default = "hsapiens",
              metavar = "ID", help = "gprofiler organism [default: %default]"),
  make_option("--direction", type = "character", default = "no",
              metavar = "yes/no",
              help = paste(
                "split ORA by regulation direction: yes = up and down separately,",
                "no = all significant genes together [default: %default]"
              )),
  make_option("--prefix", type = "character", default = "",
              metavar = "NAME",
              help = "optional prefix for output file names [default: none]")
)

opt_parser <- OptionParser(
  usage = "Rscript gprof_GO.R --in FILE --outdir DIR [--fdr 0.05] [--log2fc 1] [--direction no] [--prefix NAME]",
  option_list = option_list,
  description = paste(
    "Run GO over-representation analysis (ORA) with gprofiler2.",
    "With --direction yes, significant genes are split into up- and down-regulated sets.",
    "With --direction no, all significant genes are tested together.",
    "All GO ontologies (BP, MF, CC) are tested against the measured-gene background.",
    "",
    "Examples:",
    "  Rscript gprof_GO.R --in edger_cell_sph_vs_BHS.csv --outdir res",
    "  Rscript gprof_GO.R --in edger.csv --outdir res --fdr 0.05 --log2fc 1",
    "  Rscript gprof_GO.R --in edger.csv --outdir res --direction yes",
    "  Rscript gprof_GO.R --in young.csv --outdir res --prefix young",
    sep = "\n"
  )
)

opt <- parse_args(opt_parser)

if (is.null(opt$input_file) || is.null(opt$outdir)) {
  print_help(opt_parser)
  quit(save = "no", status = 1)
}

GO_SOURCES <- c("GO:BP", "GO:MF", "GO:CC")
ORA_FDR <- 0.05

input_file <- opt$input_file
RESULTS_DIR <- opt$outdir
fdr_cutoff <- opt$fdr_cutoff
log2fc_cutoff <- opt$log2fc_cutoff
organism <- opt$organism
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

suppressPackageStartupMessages({
  library(gprofiler2)
  library(readr)
  library(dplyr)
})

# Collapse list-columns so the enrichment table can be written as CSV.
flatten_list_cols <- function(df) {
  list_cols <- names(df)[vapply(df, is.list, logical(1))]
  for (col in list_cols) {
    df[[col]] <- vapply(
      df[[col]],
      function(x) paste(x, collapse = ","),
      character(1)
    )
  }
  df
}

# Map comma-separated Ensembl IDs to gene symbols using the DE table.
intersection_to_symbols <- function(intersection, gene_map) {
  vapply(strsplit(as.character(intersection), ",", fixed = TRUE), function(ids) {
    ids <- ids[nzchar(ids)]
    if (length(ids) == 0) {
      return("")
    }
    syms <- gene_map$gene_symbol[match(ids, gene_map$ensembl_gene)]
    paste(ifelse(is.na(syms) | syms == "", ids, syms), collapse = ",")
  }, character(1))
}

run_gost <- function(gene_ids, background, direction) {
  if (length(gene_ids) == 0) {
    message("Skipping ", direction, ": no significant genes at the chosen cutoffs.")
    return(NULL)
  }

  message(
    "Running gprofiler ORA for ", direction, " genes (n = ", length(gene_ids), ")"
  )

  gost(
    query = gene_ids,
    organism = organism,
    sources = GO_SOURCES,
    correction_method = "fdr",
    user_threshold = ORA_FDR,
    custom_bg = background,
    domain_scope = "custom_annotated",
    significant = FALSE,
    evcodes = TRUE,
    ordered_query = FALSE
  )
}

n_gost_terms <- function(gost_res, significant_only = TRUE) {
  if (is.null(gost_res) || is.null(gost_res$result) || nrow(gost_res$result) == 0) {
    return(0)
  }
  res <- gost_res$result
  if (significant_only && "significant" %in% colnames(res)) {
    return(sum(res$significant %in% TRUE))
  }
  nrow(res)
}

empty_ora_table <- function() {
  tibble(
    direction = character(),
    source = character(),
    term_id = character(),
    term_name = character(),
    p_value = double()
  )
}

significant_ora_table <- function(df) {
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }
  if ("significant" %in% colnames(df)) {
    return(df %>% filter(.data$significant %in% TRUE))
  }
  df %>% filter(p_value < ORA_FDR)
}

gost_to_table <- function(gost_res, direction, gene_map) {
  if (is.null(gost_res) || is.null(gost_res$result) || nrow(gost_res$result) == 0) {
    return(NULL)
  }

  flatten_list_cols(gost_res$result) %>%
    mutate(
      direction = direction,
      intersection_genes = intersection_to_symbols(intersection, gene_map)
    ) %>%
    relocate(direction, source, term_id, term_name, p_value)
}

# The input file
message("Input:     ", input_file)
message("Output:    ", RESULTS_DIR)
message("Prefix:    ", if (nzchar(prefix)) prefix else "(none)")
message("Organism:  ", organism)
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

background <- gene_map$ensembl_gene

deg_sig <- gene_map %>%
  filter(
    FDR < fdr_cutoff,
    abs(log2FoldChange) >= log2fc_cutoff
  )

message("Background genes:     ", length(background))
message("Total significant DEGs: ", nrow(deg_sig))

if (nrow(deg_sig) == 0) {
  stop(
    "No significant genes at FDR < ", fdr_cutoff,
    " and |log2FoldChange| >= ", log2fc_cutoff,
    call. = FALSE
  )
}

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

gene_cols <- c("ensembl_gene", "gene_symbol", "log2FoldChange", "PValue", "FDR")
gene_cols <- gene_cols[gene_cols %in% colnames(gene_map)]

ora_csv <- file.path(RESULTS_DIR, out_name("gprofiler_GO.csv"))
ora_all_csv <- file.path(RESULTS_DIR, out_name("gprofiler_GO_all.csv"))
ora_rds <- file.path(RESULTS_DIR, out_name("gprofiler_GO.rds"))

params <- list(
  input_file = input_file,
  organism = organism,
  fdr_cutoff = fdr_cutoff,
  log2fc_cutoff = log2fc_cutoff,
  direction = split_direction,
  prefix = prefix,
  ora_fdr = ORA_FDR,
  sources = GO_SOURCES,
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

  gost_up <- run_gost(up_genes$ensembl_gene, background, "up")
  gost_down <- run_gost(down_genes$ensembl_gene, background, "down")

  result_table <- bind_rows(
    gost_to_table(gost_up, "up", gene_map),
    gost_to_table(gost_down, "down", gene_map)
  )

  gost_res <- list(up = gost_up, down = gost_down, params = params)

  term_counts <- c(up = n_gost_terms(gost_up), down = n_gost_terms(gost_down))
} else {
  all_genes_csv <- file.path(RESULTS_DIR, out_name("ora_genes_all.csv"))
  write_csv(select(deg_sig, all_of(gene_cols)), all_genes_csv)
  gene_csvs <- all_genes_csv

  gost_all <- run_gost(deg_sig$ensembl_gene, background, "none")

  result_table <- gost_to_table(gost_all, "none", gene_map)
  gost_res <- list(none = gost_all, params = params)

  term_counts <- c(none = n_gost_terms(gost_all))
}

if (is.null(result_table) || nrow(result_table) == 0) {
  result_table <- empty_ora_table()
}

sig_table <- significant_ora_table(result_table)
if (nrow(sig_table) == 0) {
  if (split_direction == "yes") {
    message("No significant GO terms found for either direction.")
  } else {
    message("No significant GO terms found.")
  }
}

write_csv(sig_table, ora_csv)
write_csv(result_table, ora_all_csv)
saveRDS(gost_res, ora_rds)

message("================")
message("GO terms (gprofiler FDR < ", ORA_FDR, ")")
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
