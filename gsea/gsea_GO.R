# Usage:
#   Rscript gsea_GO.R --in edger.csv --outdir res [--ont all] [--rank signed_logp] [--prefix NAME]
#
# If --rank is omitted, create_gene_list() uses its default (signed_logp).
# --ont all runs BP, CC, and MF.
# --prefix is prepended to output file names; default is none.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--in", type = "character", default = "edger.csv", dest = "input_file",
              metavar = "FILE", help = "edgeR differential expression CSV"),
  make_option("--outdir", type = "character", default = "gsea_go",
              metavar = "DIR", help = "output directory for tables and RDS objects"),
  make_option("--ont", type = "character", default = "all",
              metavar = "BP|CC|MF|all",
              help = "GO ontology, or all for BP+CC+MF [default: %default]"),
  make_option("--organism", type = "character", default = "hsapiens",
              metavar = "ID", help = "organism: hsapiens or mmusculus [default: %default]"),
  make_option("--rank", type = "character", default = "signed_logp", dest = "ranking_method",
              metavar = "signed_logp|log2FC",
              help = "ranking statistic [default: signed_logp]"),
  make_option("--prefix", type = "character", default = "",
              metavar = "NAME",
              help = "optional prefix for output file names [default: none]")
)

opt_parser <- OptionParser(
  usage = "Rscript gsea_GO.R --in FILE --outdir DIR [--ont all] [--rank METHOD] [--prefix NAME]",
  option_list = option_list,
  description = paste(
    "Run GO GSEA with fgsea and MSigDB C5 gene sets.",
    "No plots are written; results are saved as CSV and RDS.",
    "",
    "Examples:",
    "  Rscript gsea_GO.R --in edger.csv --outdir res",
    "  Rscript gsea_GO.R --in edger.csv --outdir res --ont BP --rank log2FC",
    "  Rscript gsea_GO.R --in young.csv --outdir res --prefix young",
    sep = "\n"
  )
)

opt <- parse_args(opt_parser)

if (is.null(opt$input_file) || is.null(opt$outdir)) {
  print_help(opt_parser)
  quit(save = "no", status = 1)
}

GO_ONTS <- c("BP", "CC", "MF")
MIN_GS_SIZE <- 15
MAX_GS_SIZE <- 500
pval_cutoff <- 0.05
seed_value <- 123

input_file <- opt$input_file
RESULTS_DIR <- opt$outdir
ONT <- toupper(opt$ont)
organism <- tolower(opt$organism)
ranking_method <- opt$ranking_method
prefix <- if (is.null(opt$prefix)) "" else as.character(opt$prefix)

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

collapse_leading_edge <- function(df) {
  df %>%
    mutate(
      leadingEdge = vapply(
        leadingEdge,
        function(x) paste(x, collapse = "/"),
        character(1)
      )
    )
}

run_fgsea_ont <- function(ont, geneList, ens2sym_map, msig_org) {
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
    filter(!is.na(ensembl_gene), ensembl_gene != "", !is.na(gs_name), gs_name != "")

  if (nrow(msig) == 0) {
    stop("No MSigDB ", subcollection, " gene sets with Ensembl IDs were returned.", call. = FALSE)
  }

  gs_meta <- msig %>%
    distinct(gs_name, .keep_all = TRUE) %>%
    transmute(
      pathway = gs_name,
      ID = gs_exact_source,
      Description = str_to_sentence(gsub("_", " ", sub("^GO(BP|CC|MF)_", "", gs_name))),
      gs_id = gs_id
    )

  pathways <- lapply(split(msig$ensembl_gene, msig$gs_name), unique)
  n_in_rank <- sum(unique(unlist(pathways, use.names = FALSE)) %in% names(geneList))

  message("  gene sets: ", length(pathways))
  message("  genes overlapping the ranked list: ", n_in_rank)

  set.seed(seed_value)
  fgsea_res <- fgsea(
    pathways = pathways,
    stats = geneList,
    minSize = MIN_GS_SIZE,
    maxSize = MAX_GS_SIZE,
    nproc = 1
  )

  map_leading_edge <- function(ids) {
    if (length(ids) == 0) {
      return(character())
    }
    syms <- unname(ens2sym_map[ids])
    ifelse(is.na(syms) | syms == "", ids, syms)
  }
  fgsea_res$leadingEdge <- lapply(fgsea_res$leadingEdge, map_leading_edge)

  result_table <- as.data.frame(fgsea_res, stringsAsFactors = FALSE) %>%
    left_join(gs_meta, by = "pathway") %>%
    mutate(ONTOLOGY = ont) %>%
    relocate(ONTOLOGY, ID, Description, pathway) %>%
    arrange(padj, desc(abs(NES)))

  list(
    result = fgsea_res,
    table = result_table,
    pathways = pathways,
    db_version = unique(msig$db_version),
    subcollection = subcollection
  )
}

if (!ONT %in% c(GO_ONTS, "ALL")) {
  stop("--ont must be BP, CC, MF, or all", call. = FALSE)
}

if (!is.null(ranking_method) && !ranking_method %in% c("signed_logp", "log2FC")) {
  stop("--rank must be signed_logp or log2FC", call. = FALSE)
}

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

onts <- if (ONT == "ALL") GO_ONTS else ONT
msig_org <- organism_to_msig(organism)

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(readr)
  library(dplyr)
  library(stringr)
})

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

message("Input:     ", input_file)
message("Output:    ", RESULTS_DIR)
message("Prefix:    ", if (nzchar(prefix)) prefix else "(none)")
message("Organism:  ", organism, " (", msig_org$species, ", MSigDB ", msig_org$db_species, ")")
message("ONT:       ", paste(onts, collapse = ", "))
message("Rank:      ", ranking_method)

de <- read_csv(input_file, show_col_types = FALSE)

if (is.null(ranking_method)) {
  ranked <- create_gene_list(de)
} else {
  ranked <- create_gene_list(de, ranking_method = ranking_method)
}
geneList <- ranked$gene_list

message("Ranked genes: ", length(geneList))

ens2sym <- ranked$rank_df %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct(ensembl_gene, .keep_all = TRUE)
ens2sym_map <- setNames(ens2sym$gene_symbol, ens2sym$ensembl_gene)

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

ranked_csv <- file.path(RESULTS_DIR, out_name("ranked_genes_GO.csv"))
ranked_rds <- file.path(RESULTS_DIR, out_name("ranked_genes_GO.rds"))
write_csv(ranked$rank_df, ranked_csv)
saveRDS(geneList, ranked_rds)

ont_res <- lapply(onts, run_fgsea_ont, geneList = geneList, ens2sym_map = ens2sym_map, msig_org = msig_org)
names(ont_res) <- onts

csv_files <- character()
rds_files <- character()
sig_counts <- integer()

for (ont in onts) {
  res <- ont_res[[ont]]
  gsea_csv <- file.path(RESULTS_DIR, out_name(paste0("fgsea_GO_", ont, ".csv")))
  gsea_rds <- file.path(RESULTS_DIR, out_name(paste0("fgsea_GO_", ont, ".rds")))

  write_csv(collapse_leading_edge(res$table), gsea_csv)
  saveRDS(
    list(
      result = res$result,
      table = res$table,
      pathways = res$pathways,
      stats = geneList,
      params = list(
        input_file = input_file,
        organism = organism,
        species = msig_org$species,
        db_species = msig_org$db_species,
        ontology = ont,
        collection = "C5",
        subcollection = res$subcollection,
        ranking_method = ranking_method,
        prefix = prefix,
        minGSSize = MIN_GS_SIZE,
        maxGSSize = MAX_GS_SIZE,
        pval_cutoff = pval_cutoff,
        seed = seed_value,
        db_version = res$db_version
      )
    ),
    gsea_rds
  )

  csv_files <- c(csv_files, gsea_csv)
  rds_files <- c(rds_files, gsea_rds)
  sig_counts[[ont]] <- sum(res$table$padj < pval_cutoff, na.rm = TRUE)
}

if (length(onts) > 1) {
  combined_table <- bind_rows(lapply(ont_res, `[[`, "table")) %>%
    arrange(padj, desc(abs(NES)))
  combined_csv <- file.path(RESULTS_DIR, out_name("fgsea_GO_all.csv"))
  combined_rds <- file.path(RESULTS_DIR, out_name("fgsea_GO_all.rds"))
  write_csv(collapse_leading_edge(combined_table), combined_csv)
  saveRDS(
    list(
      by_ontology = ont_res,
      table = combined_table,
      stats = geneList,
      params = list(
        input_file = input_file,
        organism = organism,
        species = msig_org$species,
        db_species = msig_org$db_species,
        ontology = onts,
        collection = "C5",
        ranking_method = ranking_method,
        prefix = prefix,
        minGSSize = MIN_GS_SIZE,
        maxGSSize = MAX_GS_SIZE,
        pval_cutoff = pval_cutoff,
        seed = seed_value
      )
    ),
    combined_rds
  )
  csv_files <- c(csv_files, combined_csv)
  rds_files <- c(rds_files, combined_rds)
}

message("================")
message("Significant terms (padj < ", pval_cutoff, "):")
for (ont in onts) {
  n_tested <- nrow(ont_res[[ont]]$table)
  message("  ", ont, ": ", sig_counts[[ont]], " / ", n_tested, " tested")
}
message("CSV files written:")
message("  ", ranked_csv)
for (f in csv_files) {
  message("  ", f)
}
message("RDS files written:")
message("  ", ranked_rds)
for (f in rds_files) {
  message("  ", f)
}
