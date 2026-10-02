#!/usr/bin/env Rscript

usage_text <- function() {
  paste(
    "enrichviz",
    "",
    "Usage:",
    "  enrichviz <command> [script options...]",
    "  enrichviz --help",
    "",
    "Run gene-set enrichment and plots. Use -h after a command for that script's options.",
    "",
    "ORA:",
    "  enrichviz ora_gprofiler        GO over-representation with gprofiler2",
    "  enrichviz ora_barplot          top-term barplot",
    "  enrichviz ora_lollipop         top-term lollipop",
    "  enrichviz ora_dotplot          gene-ratio dotplot",
    "  enrichviz ora_upset            gene-overlap UpSet plot",
    "  enrichviz ora_ssplot           semantic-space plot (query-gene Jaccard)",
    "  enrichviz ora_simplify         GO semantic-similarity heatmap",
    "",
    "GSEA:",
    "  enrichviz gsea_go              GO GSEA with fgsea",
    "  enrichviz gsea_kegg            KEGG GSEA with fgsea",
    "  enrichviz gsea_barplot         NES barplot",
    "  enrichviz gsea_lollipop        NES lollipop",
    "  enrichviz gsea_dotplot         leading-edge-ratio dotplot",
    "  enrichviz gsea_ridgeplot       leading-edge rank-score ridges",
    "  enrichviz gsea_volcano         NES vs adjusted P",
    "  enrichviz gsea_nes             running ES for one term",
    "  enrichviz gsea_nes-compare     running ES across GSEA runs",
    "",
    "Examples:",
    "  enrichviz ora_barplot --in ora_out/gprofiler_GO.csv --outdir ora_out -n 10 --ont BP",
    "  enrichviz gsea_go --in edger.csv --outdir gsea_out --ont all",
    "",
    "From pixi without a shell, put \"--\" before the command.",
    "  pixi run enrichviz -- ora_barplot -h",
    sep = "\n"
  )
}

die <- function(...) {
  cat(paste0(...), "\n", file = stderr())
  quit(save = "no", status = 2)
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) {
  die("Could not determine the enrichviz script path.")
}
this_file <- normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE)
root <- dirname(dirname(this_file))

commands <- list(
  ora_gprofiler = list(script = "src/ora/gprof_GO.R", fn = "ora_gprofiler"),
  ora_gprof = list(script = "src/ora/gprof_GO.R", fn = "ora_gprofiler"),
  ora_barplot = list(script = "src/ora/ora_barplot.R", fn = "ora_barplot"),
  ora_lollipop = list(script = "src/ora/ora_lollipop.R", fn = "ora_lollipop"),
  ora_dotplot = list(script = "src/ora/ora_dotplot.R", fn = "ora_dotplot"),
  ora_upset = list(script = "src/ora/ora_upsetplot.R", fn = "ora_upset"),
  ora_upsetplot = list(script = "src/ora/ora_upsetplot.R", fn = "ora_upset"),
  ora_ssplot = list(script = "src/ora/ora_ssplot.R", fn = "ora_ssplot"),
  ora_simplify = list(script = "src/ora/simplify_GO.R", fn = "ora_simplify"),
  ora_simplifygo = list(script = "src/ora/simplify_GO.R", fn = "ora_simplify"),
  "ora_simplify-go" = list(script = "src/ora/simplify_GO.R", fn = "ora_simplify"),
  gsea_go = list(script = "src/gsea/gsea_GO.R", fn = "gsea_go"),
  gsea_kegg = list(script = "src/gsea/gsea_KEGG.R", fn = "gsea_kegg"),
  gsea_barplot = list(script = "src/gsea/gsea_barplot.R", fn = "gsea_barplot"),
  gsea_lollipop = list(script = "src/gsea/gsea_lollipop.R", fn = "gsea_lollipop"),
  gsea_dotplot = list(script = "src/gsea/gsea_dotplot.R", fn = "gsea_dotplot"),
  gsea_ridgeplot = list(script = "src/gsea/gsea_ridgeplot.R", fn = "gsea_ridgeplot"),
  gsea_volcano = list(script = "src/gsea/gsea_volcano.R", fn = "gsea_volcano"),
  gsea_nes = list(script = "src/gsea/gsea_NES.R", fn = "gsea_nes"),
  "gsea_nes-compare" = list(script = "src/gsea/gsea_NES_compare.R", fn = "gsea_nes_compare"),
  gsea_nes_compare = list(script = "src/gsea/gsea_NES_compare.R", fn = "gsea_nes_compare")
)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0L) {
  cat(usage_text(), "\n", sep = "")
  quit(save = "no", status = 2)
}

cmd <- args[[1]]
cmd_args <- args[-1]

if (cmd %in% c("-h", "--help", "help")) {
  cat(usage_text(), "\n", sep = "")
  quit(save = "no", status = 0)
}

if (cmd %in% c("ora", "gsea")) {
  die(
    "Commands are a single word, for example: enrichviz ", cmd, "_barplot\n",
    usage_text()
  )
}

spec <- commands[[cmd]]

if (is.null(spec)) {
  die("Unknown command: ", cmd, "\n", usage_text())
}

script <- file.path(root, spec$script)
if (!file.exists(script)) {
  die("Script not found: ", script)
}

source(file.path(root, "src", "cli.R"))
ENRICHVIZ_MODULE <- normalizePath(script, winslash = "/", mustWork = TRUE)
source(script, local = FALSE)

options_fn <- get0(paste0(spec$fn, "_options"), mode = "function", inherits = TRUE)
run_fn <- get0(paste0(spec$fn, "_run"), mode = "function", inherits = TRUE)
if (is.null(options_fn) || is.null(run_fn)) {
  die("Module is missing ", spec$fn, "_options() or ", spec$fn, "_run()")
}
invoke_module(options_fn, run_fn, args = cmd_args)
