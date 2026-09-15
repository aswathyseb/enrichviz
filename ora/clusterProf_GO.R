# Usage:
#   Rscript clusterProf_GO.R --in edger.csv --outdir res [--fdr 0.05] [--log2fc 1] [--direction no]
#
# Select significant DEGs by FDR and |log2FC|, then run clusterProfiler
# enrichGO across GO BP, MF, and CC. With --direction yes, up- and
# down-regulated genes are tested separately; with --direction no
# (default), all significant genes are tested together.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--in", type = "character", default = "edger.csv", dest = "input_file",
              metavar = "FILE", help = "edgeR differential expression CSV"),
  make_option("--outdir", type = "character", default = "clusterProf_go",
              metavar = "DIR", help = "output directory for tables and RDS objects"),
  make_option("--fdr", type = "double", default = 0.05, dest = "fdr_cutoff",
              metavar = "NUM", help = "DEG FDR cutoff [default: %default]"),
  make_option("--log2fc", type = "double", default = 1, dest = "log2fc_cutoff",
              metavar = "NUM", help = "DEG |log2FoldChange| cutoff [default: %default]"),
  make_option("--organism", type = "character", default = "hsapiens",
              metavar = "ID", help = "organism: hsapiens or mmusculus [default: %default]"),
  make_option("--direction", type = "character", default = "no",
              metavar = "yes/no",
              help = paste(
                "split ORA by regulation direction: yes = up and down separately,",
                "no = all significant genes together [default: %default]"
              ))
)

opt_parser <- OptionParser(
  usage = "Rscript clusterProf_GO.R --in FILE --outdir DIR [--fdr 0.05] [--log2fc 1] [--direction no]",
  option_list = option_list,
  description = paste(
    "Run GO over-representation analysis (ORA) with clusterProfiler::enrichGO.",
    "With --direction yes, significant genes are split into up- and down-regulated sets.",
    "With --direction no, all significant genes are tested together.",
    "All GO ontologies (BP, MF, CC) are tested against the measured-gene background.",
    "",
    "Examples:",
    "  Rscript clusterProf_GO.R --in edger_cell_sph_vs_BHS.csv --outdir res",
    "  Rscript clusterProf_GO.R --in edger.csv --outdir res --fdr 0.05 --log2fc 1",
    "  Rscript clusterProf_GO.R --in edger.csv --outdir res --direction yes",
    sep = "\n"
  )
)

opt <- parse_args(opt_parser)

if (is.null(opt$input_file) || is.null(opt$outdir)) {
  print_help(opt_parser)
  quit(save = "no", status = 1)
}

GO_ONTS <- c("BP", "MF", "CC")
ORA_FDR <- 0.05

input_file <- opt$input_file
RESULTS_DIR <- opt$outdir
fdr_cutoff <- opt$fdr_cutoff
log2fc_cutoff <- opt$log2fc_cutoff
organism <- tolower(opt$organism)
split_direction <- tolower(opt$direction)

organism_to_orgdb <- function(organism) {
  orgdb_map <- c(
    hsapiens = "org.Hs.eg.db",
    human = "org.Hs.eg.db",
    hs = "org.Hs.eg.db",
    mmusculus = "org.Mm.eg.db",
    mouse = "org.Mm.eg.db",
    mm = "org.Mm.eg.db"
  )
  orgdb <- unname(orgdb_map[organism])
  if (length(orgdb) != 1 || is.na(orgdb)) {
    stop("--organism must be hsapiens or mmusculus", call. = FALSE)
  }
  orgdb
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

orgdb_name <- organism_to_orgdb(organism)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
  library(readr)
  library(dplyr)
  library(tibble)
})

if (!requireNamespace(orgdb_name, quietly = TRUE)) {
  stop("Install the OrgDb package ", orgdb_name, call. = FALSE)
}
suppressPackageStartupMessages(library(orgdb_name, character.only = TRUE))
OrgDb <- get(orgdb_name)

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

n_enrich_terms <- function(enrich_list, significant_only = TRUE) {
  if (is.null(enrich_list)) {
    return(0)
  }
  sum(vapply(enrich_list, function(ego) {
    if (is.null(ego)) {
      return(0)
    }
    df <- as.data.frame(ego)
    if (nrow(df) == 0) {
      return(0)
    }
    if (significant_only) {
      return(sum(df$p.adjust < ORA_FDR, na.rm = TRUE))
    }
    nrow(df)
  }, numeric(1)))
}

enrich_to_table <- function(enrich_list, direction) {
  if (is.null(enrich_list)) {
    return(NULL)
  }

  pieces <- lapply(names(enrich_list), function(ont) {
    ego <- enrich_list[[ont]]
    if (is.null(ego)) {
      return(NULL)
    }
    df <- as.data.frame(ego)
    if (nrow(df) == 0) {
      return(NULL)
    }
    if (!"ONTOLOGY" %in% colnames(df)) {
      df$ONTOLOGY <- ont
    }
    df$direction <- direction
    df
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (length(pieces) == 0) {
    return(NULL)
  }

  bind_rows(pieces) %>%
    relocate(direction, ONTOLOGY, ID, Description, p.adjust)
}

run_enrichGO <- function(gene_ids, background, direction) {
  gene_ids <- unique(gene_ids)
  background <- unique(background)

  if (length(gene_ids) == 0) {
    message("Skipping ", direction, ": no significant genes at the chosen cutoffs.")
    return(NULL)
  }

  message(
    "Running clusterProfiler ORA for ", direction,
    " genes (n = ", length(gene_ids), ")"
  )

  res <- lapply(GO_ONTS, function(ont) {
    message("  ontology ", ont)
    tryCatch(
      enrichGO(
        gene = gene_ids,
        universe = background,
        OrgDb = OrgDb,
        keyType = "ENSEMBL",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        minGSSize = 10,
        maxGSSize = Inf,
        readable = TRUE
      ),
      error = function(e) {
        message("  enrichGO failed for ", ont, ": ", conditionMessage(e))
        NULL
      }
    )
  })
  names(res) <- GO_ONTS
  res
}

message("Input:     ", input_file)
message("Output:    ", RESULTS_DIR)
message("Organism:  ", organism, " (", orgdb_name, ")")
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

ora_csv <- file.path(RESULTS_DIR, "clusterProf_GO.csv")
ora_all_csv <- file.path(RESULTS_DIR, "clusterProf_GO_all.csv")
ora_rds <- file.path(RESULTS_DIR, "clusterProf_GO.rds")

params <- list(
  input_file = input_file,
  organism = organism,
  orgdb = orgdb_name,
  fdr_cutoff = fdr_cutoff,
  log2fc_cutoff = log2fc_cutoff,
  direction = split_direction,
  ora_fdr = ORA_FDR,
  ontologies = GO_ONTS,
  keyType = "ENSEMBL",
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

  up_genes_csv <- file.path(RESULTS_DIR, "ora_genes_up.csv")
  down_genes_csv <- file.path(RESULTS_DIR, "ora_genes_down.csv")
  write_csv(select(up_genes, all_of(gene_cols)), up_genes_csv)
  write_csv(select(down_genes, all_of(gene_cols)), down_genes_csv)
  gene_csvs <- c(up_genes_csv, down_genes_csv)

  enrich_up <- run_enrichGO(up_genes$ensembl_gene, background, "up")
  enrich_down <- run_enrichGO(down_genes$ensembl_gene, background, "down")

  result_table <- bind_rows(
    enrich_to_table(enrich_up, "up"),
    enrich_to_table(enrich_down, "down")
  )

  enrich_res <- list(up = enrich_up, down = enrich_down, params = params)
  term_counts <- c(up = n_enrich_terms(enrich_up), down = n_enrich_terms(enrich_down))
} else {
  all_genes_csv <- file.path(RESULTS_DIR, "ora_genes_all.csv")
  write_csv(select(deg_sig, all_of(gene_cols)), all_genes_csv)
  gene_csvs <- all_genes_csv

  enrich_all <- run_enrichGO(deg_sig$ensembl_gene, background, "none")

  result_table <- enrich_to_table(enrich_all, "none")
  enrich_res <- list(none = enrich_all, params = params)
  term_counts <- c(none = n_enrich_terms(enrich_all))
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
saveRDS(enrich_res, ora_rds)

message("================")
message("GO terms (clusterProfiler BH FDR < ", ORA_FDR, ")")
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
