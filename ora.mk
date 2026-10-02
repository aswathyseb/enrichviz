# ORA plots from a differential-expression table or an existing ORA result.
#
#   make -f ora.mk
#   make -f ora.mk IN=data/edger.csv OUT=res ONT=BP DIRECTION=up
#   make -f ora.mk ORA=res/gprofiler_GO.csv OUT=res ONT=BP
#
# IN is a differential-expression table. ora_gprofiler runs, then the plots.
# ORA is an existing enrichment table. The plots use it and gprofiler is skipped.

IN        ?= data/edger.csv
OUT       ?= res
ONT       ?= BP
N         ?= 10
N_SS      ?= 30
DIRECTION ?= up
ORGANISM  ?= hsapiens
FDR       ?= 0.05
LOG2FC    ?= 1
ORA       ?=

ENRICHVIZ := Rscript src/enrichviz.R

# ORA names an enrichment table. Otherwise IN is the differential-expression file.
ifneq (${ORA},)
  ORA_TABLE := ${ORA}
  FROM_DE :=
  ifeq ($(wildcard ${ORA_TABLE}),)
    $(error ORA result not found: ${ORA_TABLE})
  endif
else
  ORA_TABLE := ${OUT}/gprofiler_GO.csv
  FROM_DE := yes
endif

# Create separate plots for up and down directions
ifeq (${DIRECTION},none)
  UPSET_PDF := ${OUT}/ora_${ONT}_upset.pdf
  SS_PDF    := ${OUT}/ora_${ONT}_ssplot.pdf
else
  UPSET_PDF := ${OUT}/ora_${ONT}_${DIRECTION}_upset.pdf
  SS_PDF    := ${OUT}/ora_${ONT}_${DIRECTION}_ssplot.pdf
endif

# Up and Down are added to the same plot when direction is specified
BAR_PDF      := ${OUT}/ora_${ONT}_barplot.pdf
LOLLIPOP_PDF := ${OUT}/ora_${ONT}_lollipop.pdf
DOT_PDF      := ${OUT}/ora_${ONT}_dotplot.pdf
SIMPLIFY_PDF := ${OUT}/GO_${ONT}_${DIRECTION}_simplifyGO.pdf

PLOTS := ${BAR_PDF} ${LOLLIPOP_PDF} ${DOT_PDF} ${UPSET_PDF} ${SS_PDF} ${SIMPLIFY_PDF}

.PHONY: all plots gprofiler barplot lollipopdotplot upsetplot ssplot simplifygo_plot clean

all: plots

plots: ${PLOTS}
	ls -l ${PLOTS}

gprofiler: ${ORA_TABLE}

barplot: ${BAR_PDF}
	ls -l ${BAR_PDF}

lollipop: ${LOLLIPOP_PDF}
	ls -l ${LOLLIPOP_PDF}

dotplot: ${DOT_PDF}
	ls -l ${DOT_PDF}

upsetplot: ${UPSET_PDF}
	ls -l ${UPSET_PDF}

ssplot: ${SS_PDF}
	ls -l ${SS_PDF}

simplifygo_plots: ${SIMPLIFY_PDF}
	ls -l ${SIMPLIFY_PDF}

ifeq (${FROM_DE},yes)
${ORA_TABLE}: ${IN}
	mkdir -p ${OUT}
	${ENRICHVIZ} ora_gprofiler \
	  --in ${IN} --outdir ${OUT} \
	  --direction yes --organism ${ORGANISM} \
	  --fdr ${FDR} --log2fc ${LOG2FC}
endif

${BAR_PDF}: ${ORA_TABLE}
	${ENRICHVIZ} ora_barplot --in ${ORA_TABLE} --outdir ${OUT} -n ${N} --ont ${ONT}

${LOLLIPOP_PDF}: ${ORA_TABLE}
	${ENRICHVIZ} ora_lollipop --in ${ORA_TABLE} --outdir ${OUT} -n ${N} --ont ${ONT}

${DOT_PDF}: ${ORA_TABLE}
	${ENRICHVIZ} ora_dotplot --in ${ORA_TABLE} --outdir ${OUT} -n ${N} --ont ${ONT}

${UPSET_PDF}: ${ORA_TABLE}
	${ENRICHVIZ} ora_upset --in ${ORA_TABLE} --outdir ${OUT} -n ${N} --ont ${ONT} --direction ${DIRECTION}

${SS_PDF}: ${ORA_TABLE}
	${ENRICHVIZ} ora_ssplot --in ${ORA_TABLE} --outdir ${OUT} -n ${N_SS} --ont ${ONT} --direction ${DIRECTION}

${SIMPLIFY_PDF}: ${ORA_TABLE}
	${ENRICHVIZ} ora_simplify --in ${ORA_TABLE} --outdir ${OUT} --ont ${ONT} --direction ${DIRECTION} --organism ${ORGANISM}

clean:
	rm -f \
	  $(if ${FROM_DE},${OUT}/gprofiler_GO.csv ${OUT}/gprofiler_GO_all.csv ${OUT}/gprofiler_GO.rds ${OUT}/ora_genes_up.csv ${OUT}/ora_genes_down.csv,) \
	  ${BAR_PDF} ${OUT}/ora_${ONT}_barplot_table.csv \
	  ${LOLLIPOP_PDF} ${OUT}/ora_${ONT}_lollipop_table.csv \
	  ${DOT_PDF} ${OUT}/ora_${ONT}_dotplot_table.csv \
	  ${UPSET_PDF} ${UPSET_PDF:.pdf=_table.csv} \
	  ${SS_PDF} ${SS_PDF:.pdf=_table.csv} ${SS_PDF:.pdf=_similarity.rds} \
	  ${SIMPLIFY_PDF} \
	  ${OUT}/GO_${ONT}_${DIRECTION}_simplifyGO_clusters.csv \
	  ${OUT}/GO_${ONT}_${DIRECTION}_similarity.rds
