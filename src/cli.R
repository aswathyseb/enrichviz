# Shared command-line options for enrichviz scripts.
# Each script sources this file and passes only the options it uses.

suppressPackageStartupMessages(library(optparse))

.enrichviz_parser <- NULL

opt_in <- function(default, help = "input CSV/TSV or RDS") {
  make_option("--in", type = "character", default = default, dest = "input_file",
              metavar = "FILE", help = help)
}

opt_outdir <- function(default, help = "output directory") {
  make_option("--outdir", type = "character", default = default,
              metavar = "DIR", help = help)
}

opt_n <- function(default = 10L, help = "top terms [default: %default]") {
  make_option(c("-n", "--n"), type = "integer", default = as.integer(default),
              dest = "top_n", metavar = "N", help = help)
}

opt_ont <- function(default = "BP",
                    metavar = "BP|CC|MF|KEGG|all",
                    help = paste(
                      "BP, CC, or MF plots that GO ontology.",
                      "KEGG plots KEGG ids. all plots every GO id",
                      "[default: %default]"
                    )) {
  make_option("--ont", type = "character", default = default,
              metavar = metavar, help = help)
}

opt_prefix <- function(default = "") {
  make_option("--prefix", type = "character", default = default,
              metavar = "NAME",
              help = "optional prefix for output file names [default: none]")
}

opt_organism <- function(default = "hsapiens",
                         metavar = "ID",
                         help = "organism: hsapiens or mmusculus [default: %default]") {
  make_option("--organism", type = "character", default = default,
              metavar = metavar, help = help)
}

opt_fdr <- function(default = 0.05) {
  make_option("--fdr", type = "double", default = default, dest = "fdr_cutoff",
              metavar = "NUM", help = "DEG FDR cutoff [default: %default]")
}

opt_log2fc <- function(default = 1) {
  make_option("--log2fc", type = "double", default = default, dest = "log2fc_cutoff",
              metavar = "NUM", help = "DEG |log2FoldChange| cutoff [default: %default]")
}

opt_direction <- function(default, metavar = "DIR", help) {
  make_option("--direction", type = "character", default = default,
              metavar = metavar, help = help)
}

opt_rank <- function(default,
                     dest = "ranking_method",
                     metavar = "METHOD",
                     help = "ranking statistic [default: signed_logp]") {
  make_option("--rank", type = "character", default = default, dest = dest,
              metavar = metavar, help = help)
}

opt_go <- function(metavar = "ID",
                  help = "GO or KEGG term to plot (required unless --term)") {
  make_option("--GO", type = "character", default = NULL, dest = "GO_ID",
              metavar = metavar, help = help)
}

opt_term <- function(help = "alias for --GO") {
  make_option("--term", type = "character", default = NULL, dest = "term_id",
              metavar = "ID", help = help)
}

# Quit with the help text from the most recent parse_cli() call.
cli_usage <- function() {
  if (is.null(.enrichviz_parser)) {
    stop("No command-line parser has been created", call. = FALSE)
  }
  print_help(.enrichviz_parser)
  quit(save = "no", status = 1)
}

module_dir <- function() {
  mod <- get0("ENRICHVIZ_MODULE", ifnotfound = "", inherits = TRUE)
  if (is.character(mod) && length(mod) == 1L && nzchar(mod)) {
    return(dirname(normalizePath(mod, winslash = "/", mustWork = TRUE)))
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) != 1L) {
    stop("Could not determine the script directory", call. = FALSE)
  }
  dirname(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE))
}

source_sibling <- function(filename) {
  path <- file.path(module_dir(), filename)
  if (!file.exists(path)) {
    stop("Could not find ", filename, call. = FALSE)
  }
  source(path, local = FALSE)
}

invoke_module <- function(options_fn, run_fn, args = commandArgs(trailingOnly = TRUE)) {
  spec <- options_fn()
  require_input <- if (is.null(spec$require_input)) TRUE else isTRUE(spec$require_input)
  opt <- parse_cli(
    usage = spec$usage,
    description = spec$description,
    options = spec$options,
    require_input = require_input,
    args = args
  )
  run_fn(opt)
}

parse_cli <- function(usage, description, options, require_input = TRUE,
                      args = commandArgs(trailingOnly = TRUE)) {
  parser <- OptionParser(
    usage = usage,
    option_list = options,
    description = description
  )
  .enrichviz_parser <<- parser
  opt <- parse_args(parser, args = args)
  if (isTRUE(require_input) && (is.null(opt$input_file) || is.null(opt$outdir))) {
    cli_usage()
  }
  if (!is.null(opt$top_n)) {
    top_n <- suppressWarnings(as.integer(opt$top_n))
    if (length(top_n) != 1L || !is.finite(top_n) || top_n < 1L) {
      stop("-n must be a positive integer", call. = FALSE)
    }
  }
  opt
}
