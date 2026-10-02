# GSEA plots from a differential-expression table.
#
# Run from the repository root, inside a pixi shell so Rscript has the
# enrichviz packages:
#
#   make -f gsea.mk
#   make -f gsea.mk IN=data/edger.csv OUT=res ONT=BP GO_ID=GO:0051607
#
# ONT is BP, CC, or MF. It selects which GO terms are plotted.
# gsea_go still tests BP, CC, and MF so one result table feeds every ontology.
# Bar, lollipop, dot, and ridge plots show the top N terms in each NES direction.
# The volcano plot draws every tested term and labels the top N per direction.
#
# gsea_nes draws one term. Set GO_ID to a GO id, or leave it empty to use the
# first term in the GO results table. The expression table is
# passed as --expr so a leading-edge heatmap is added when the columns exist.
# gsea_nes-compare needs two GSEA runs, so it is not part of this file.

IN        ?= data/edger.csv
OUT       ?= res
ONT       ?= BP
N         ?= 10
ORGANISM  ?= hsapiens
RANK      ?= signed_logp
GO_ID     ?=

ENRICHVIZ := Rscript src/enrichviz.R

GO_CSV    := ${OUT}/fgsea_GO_all.csv
GO_RDS    := ${OUT}/fgsea_GO_all.rds
GO_RANK   := ${OUT}/ranked_genes_GO.csv
KEGG_CSV  := ${OUT}/fgsea_KEGG.csv
KEGG_RANK := ${OUT}/ranked_genes_KEGG.csv
NES_STAMP := ${OUT}/gsea_nes.stamp

BAR_PDF      := ${OUT}/gsea_${ONT}_barplot.pdf
LOLLIPOP_PDF := ${OUT}/gsea_${ONT}_lollipop.pdf
DOT_PDF      := ${OUT}/gsea_${ONT}_dotplot.pdf
RIDGE_PDF    := ${OUT}/gsea_${ONT}_ridgeplot.pdf
VOLCANO_PDF  := ${OUT}/gsea_${ONT}_volcano.pdf

KEGG_BAR_PDF      := ${OUT}/gsea_KEGG_barplot.pdf
KEGG_LOLLIPOP_PDF := ${OUT}/gsea_KEGG_lollipop.pdf
KEGG_DOT_PDF      := ${OUT}/gsea_KEGG_dotplot.pdf
KEGG_RIDGE_PDF    := ${OUT}/gsea_KEGG_ridgeplot.pdf
KEGG_VOLCANO_PDF  := ${OUT}/gsea_KEGG_volcano.pdf

GO_PLOTS := ${BAR_PDF} ${LOLLIPOP_PDF} ${DOT_PDF} ${RIDGE_PDF} ${VOLCANO_PDF}
KEGG_PLOTS := ${KEGG_BAR_PDF} ${KEGG_LOLLIPOP_PDF} ${KEGG_DOT_PDF} ${KEGG_RIDGE_PDF} ${KEGG_VOLCANO_PDF}

.PHONY: all plots gsea_go gsea_kegg \
	barplot lollipop dotplot ridgeplot volcano nes \
	kegg kegg_barplot kegg_lollipop kegg_dotplot kegg_ridgeplot kegg_volcano \
	clean

all: plots nes

plots: ${GO_PLOTS} ${KEGG_PLOTS}
	ls -l ${GO_PLOTS} ${KEGG_PLOTS}

gsea_go: ${GO_CSV}

gsea_kegg: ${KEGG_CSV}

go_plots: ${GO_PLOTS} nes
	ls -l ${GO_PLOTS}

kegg_plots: ${KEGG_PLOTS} 
	ls -l ${KEGG_PLOTS}

barplot: ${BAR_PDF}
	ls -l ${BAR_PDF}

lollipop: ${LOLLIPOP_PDF}
	ls -l ${LOLLIPOP_PDF}

dotplot: ${DOT_PDF}
	ls -l ${DOT_PDF}

ridgeplot: ${RIDGE_PDF}
	ls -l ${RIDGE_PDF}

volcano: ${VOLCANO_PDF}
	ls -l ${VOLCANO_PDF}

nes: ${NES_STAMP}

kegg: ${KEGG_PLOTS}
	ls -l ${KEGG_PLOTS}

kegg_barplot: ${KEGG_BAR_PDF}
	ls -l ${KEGG_BAR_PDF}

kegg_lollipop: ${KEGG_LOLLIPOP_PDF}
	ls -l ${KEGG_LOLLIPOP_PDF}

kegg_dotplot: ${KEGG_DOT_PDF}
	ls -l ${KEGG_DOT_PDF}

kegg_ridgeplot: ${KEGG_RIDGE_PDF}
	ls -l ${KEGG_RIDGE_PDF}

kegg_volcano: ${KEGG_VOLCANO_PDF}
	ls -l ${KEGG_VOLCANO_PDF}

${GO_CSV}: ${IN}
	mkdir -p ${OUT}
	${ENRICHVIZ} gsea_go \
	  --in ${IN} --outdir ${OUT} \
	  --ont all --organism ${ORGANISM} --rank ${RANK}

${KEGG_CSV}: ${IN}
	mkdir -p ${OUT}
	${ENRICHVIZ} gsea_kegg \
	  --in ${IN} --outdir ${OUT} \
	  --organism ${ORGANISM} --rank ${RANK}

${BAR_PDF}: ${GO_CSV}
	${ENRICHVIZ} gsea_barplot --in ${GO_CSV} --outdir ${OUT} -n ${N} --ont ${ONT}

${LOLLIPOP_PDF}: ${GO_CSV}
	${ENRICHVIZ} gsea_lollipop --in ${GO_CSV} --outdir ${OUT} -n ${N} --ont ${ONT}

${DOT_PDF}: ${GO_CSV}
	${ENRICHVIZ} gsea_dotplot --in ${GO_CSV} --outdir ${OUT} -n ${N} --ont ${ONT}

${RIDGE_PDF}: ${GO_CSV}
	${ENRICHVIZ} gsea_ridgeplot --in ${GO_CSV} --rank ${GO_RANK} --outdir ${OUT} -n ${N} --ont ${ONT}

${VOLCANO_PDF}: ${GO_CSV}
	${ENRICHVIZ} gsea_volcano --in ${GO_CSV} --outdir ${OUT} -n ${N} --ont ${ONT}

${KEGG_BAR_PDF}: ${KEGG_CSV}
	${ENRICHVIZ} gsea_barplot --in ${KEGG_CSV} --outdir ${OUT} -n ${N} --ont KEGG

${KEGG_LOLLIPOP_PDF}: ${KEGG_CSV}
	${ENRICHVIZ} gsea_lollipop --in ${KEGG_CSV} --outdir ${OUT} -n ${N} --ont KEGG

${KEGG_DOT_PDF}: ${KEGG_CSV}
	${ENRICHVIZ} gsea_dotplot --in ${KEGG_CSV} --outdir ${OUT} -n ${N} --ont KEGG

${KEGG_RIDGE_PDF}: ${KEGG_CSV}
	${ENRICHVIZ} gsea_ridgeplot --in ${KEGG_CSV} --rank ${KEGG_RANK} --outdir ${OUT} -n ${N} --ont KEGG

${KEGG_VOLCANO_PDF}: ${KEGG_CSV}
	${ENRICHVIZ} gsea_volcano --in ${KEGG_CSV} --outdir ${OUT} -n ${N} --ont KEGG

${NES_STAMP}: ${GO_CSV}
	term='${GO_ID}'; \
	if [ -z "$$term" ]; then \
	  term=$$(awk -F, 'NR == 2 { print $$2; exit }' ${GO_CSV}); \
	fi; \
	test -n "$$term"; \
	${ENRICHVIZ} gsea_nes --in ${GO_RDS} --GO "$$term" --expr ${IN} --outdir ${OUT}; \
	touch $@

clean:
	rm -f \
	  ${NES_STAMP} \
	  ${GO_CSV} ${GO_RDS} \
	  ${OUT}/fgsea_GO_BP.csv ${OUT}/fgsea_GO_BP.rds \
	  ${OUT}/fgsea_GO_CC.csv ${OUT}/fgsea_GO_CC.rds \
	  ${OUT}/fgsea_GO_MF.csv ${OUT}/fgsea_GO_MF.rds \
	  ${GO_RANK} ${OUT}/ranked_genes_GO.rds \
	  ${KEGG_CSV} ${OUT}/fgsea_KEGG.rds \
	  ${KEGG_RANK} ${OUT}/ranked_genes_KEGG.rds \
	  ${BAR_PDF} ${OUT}/gsea_${ONT}_barplot_table.csv \
	  ${LOLLIPOP_PDF} ${OUT}/gsea_${ONT}_lollipop_table.csv \
	  ${DOT_PDF} ${OUT}/gsea_${ONT}_dotplot_table.csv \
	  ${RIDGE_PDF} ${OUT}/gsea_${ONT}_ridgeplot_table.csv \
	  ${VOLCANO_PDF} ${OUT}/gsea_${ONT}_volcano_table.csv \
	  ${KEGG_BAR_PDF} ${OUT}/gsea_KEGG_barplot_table.csv \
	  ${KEGG_LOLLIPOP_PDF} ${OUT}/gsea_KEGG_lollipop_table.csv \
	  ${KEGG_DOT_PDF} ${OUT}/gsea_KEGG_dotplot_table.csv \
	  ${KEGG_RIDGE_PDF} ${OUT}/gsea_KEGG_ridgeplot_table.csv \
	  ${KEGG_VOLCANO_PDF} ${OUT}/gsea_KEGG_volcano_table.csv \
	  ${OUT}/gsea_*_NES.pdf ${OUT}/gsea_*_NES_table.csv
