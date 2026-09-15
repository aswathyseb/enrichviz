# Usage:
#   Rscript gsea_KEGG.R --in edger.csv --outdir res [--kegg legacy] [--rank signed_logp] [--prefix NAME]
#
# If --rank is omitted, create_gene_list() uses its default (signed_logp).
# --kegg legacy uses MSigDB C2 CP:KEGG_LEGACY (classic hsa##### pathway IDs).
# --kegg medicus uses MSigDB C2 CP:KEGG_MEDICUS.
# --prefix is prepended to output file names; default is none.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--in", type = "character", default = "edger.csv", dest = "input_file",
              metavar = "FILE", help = "edgeR differential expression CSV"),
  make_option("--outdir", type = "character", default = "gsea_kegg",
              metavar = "DIR", help = "output directory for tables and RDS objects"),
  make_option("--kegg", type = "character", default = "legacy", dest = "kegg_db",
              metavar = "legacy|medicus",
              help = "MSigDB KEGG collection [default: %default]"),
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
  usage = "Rscript gsea_KEGG.R --in FILE --outdir DIR [--kegg legacy] [--rank METHOD] [--prefix NAME]",
  option_list = option_list,
  description = paste(
    "Run KEGG GSEA with fgsea and MSigDB C2 gene sets.",
    "No plots are written; results are saved as CSV and RDS.",
    "KEGG collections are human MSigDB sets; mouse uses ortholog mapping.",
    "",
    "Examples:",
    "  Rscript gsea_KEGG.R --in edger.csv --outdir res",
    "  Rscript gsea_KEGG.R --in edger.csv --outdir res --kegg medicus --rank log2FC",
    "  Rscript gsea_KEGG.R --in young.csv --outdir res --prefix young",
    sep = "\n"
  )
)

opt <- parse_args(opt_parser)

if (is.null(opt$input_file) || is.null(opt$outdir)) {
  print_help(opt_parser)
  quit(save = "no", status = 1)
}

MIN_GS_SIZE <- 15
MAX_GS_SIZE <- 500
pval_cutoff <- 0.05
seed_value <- 123

kegg_collections <- list(
  legacy = list(subcollection = "CP:KEGG_LEGACY", label = "KEGG Legacy"),
  medicus = list(subcollection = "CP:KEGG_MEDICUS", label = "KEGG Medicus")
)

input_file <- opt$input_file
RESULTS_DIR <- opt$outdir
kegg_db <- tolower(opt$kegg_db)
organism <- tolower(opt$organism)
ranking_method <- opt$ranking_method
prefix <- if (is.null(opt$prefix)) "" else as.character(opt$prefix)

if (!kegg_db %in% names(kegg_collections)) {
  stop("--kegg must be legacy or medicus", call. = FALSE)
}
kegg_info <- kegg_collections[[kegg_db]]

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

kegg_description <- function(gs_name) {
  stripped <- sub("^KEGG_(MEDICUS_)?", "", gs_name)
  str_to_sentence(gsub("_", " ", stripped))
}

run_fgsea_kegg <- function(geneList, ens2sym_map, msig_org, kegg_info) {
  message("Loading MSigDB ", kegg_info$subcollection, " gene sets...")

  msig <- tryCatch(
    msigdbr(
      db_species = "HS",
      species = msig_org$species,
      collection = "C2",
      subcollection = kegg_info$subcollection
    ),
    error = function(e) {
      stop("Failed to load MSigDB gene sets: ", conditionMessage(e), call. = FALSE)
    }
  )

  msig <- msig %>%
    filter(!is.na(ensembl_gene), ensembl_gene != "", !is.na(gs_name), gs_name != "")

  if (nrow(msig) == 0) {
    stop(
      "No MSigDB ", kegg_info$subcollection,
      " gene sets with Ensembl IDs were returned.",
      call. = FALSE
    )
  }

  gs_meta <- msig %>%
    distinct(gs_name, .keep_all = TRUE) %>%
    transmute(
      pathway = gs_name,
      ID = ifelse(
        !is.na(gs_exact_source) & gs_exact_source != "",
        gs_exact_source,
        gs_id
      ),
      Description = kegg_description(gs_name),
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
    mutate(ONTOLOGY = "KEGG") %>%
    relocate(ONTOLOGY, ID, Description, pathway) %>%
    arrange(padj, desc(abs(NES)))

  list(
    result = fgsea_res,
    table = result_table,
    pathways = pathways,
    db_version = unique(msig$db_version),
    subcollection = kegg_info$subcollection
  )
}

if (!is.null(ranking_method) && !ranking_method %in% c("signed_logp", "log2FC")) {
  stop("--rank must be signed_logp or log2FC", call. = FALSE)
}

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

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
message("Organism:  ", organism, " (", msig_org$species, ")")
if (msig_org$db_species != "HS") {
  message(
    "KEGG sets: human MSigDB ", kegg_info$label,
    " mapped to ", msig_org$species, " orthologs"
  )
} else {
  message("KEGG sets: MSigDB ", kegg_info$label, " (", kegg_info$subcollection, ")")
}
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

ranked_csv <- file.path(RESULTS_DIR, gsea_out_name("ranked_genes_KEGG.csv", prefix))
ranked_rds <- file.path(RESULTS_DIR, gsea_out_name("ranked_genes_KEGG.rds", prefix))
write_csv(ranked$rank_df, ranked_csv)
saveRDS(geneList, ranked_rds)

res <- run_fgsea_kegg(geneList, ens2sym_map, msig_org, kegg_info)

gsea_csv <- file.path(RESULTS_DIR, gsea_out_name("fgsea_KEGG.csv", prefix))
gsea_rds <- file.path(RESULTS_DIR, gsea_out_name("fgsea_KEGG.rds", prefix))

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
      db_species = "HS",
      ontology = "KEGG",
      collection = "C2",
      subcollection = res$subcollection,
      kegg = kegg_db,
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

n_tested <- nrow(res$table)
n_sig <- sum(res$table$padj < pval_cutoff, na.rm = TRUE)

message("================")
message("Significant KEGG pathways (padj < ", pval_cutoff, "): ", n_sig, " / ", n_tested, " tested")
message("CSV files written:")
message("  ", ranked_csv)
message("  ", gsea_csv)
message("RDS files written:")
message("  ", ranked_rds)
message("  ", gsea_rds)
