# Input file for differentially expressed genes
IN   ?= data/edger.csv
# Output directory
OUT  ?= res
# Ontology to plot
ONT  ?= BP

# Functional analysis results
ORA_IN  ?= data/gprof.csv

PREFIX ?=


#
# Start from differential expression analysis
#

# Run gprofiler and create all ora plots
ora:
	make -f ora.mk IN=${IN} OUT=${OUT} ONT=${ONT} PREFIX=${PREFIX} plots

# Run fgsea and create all gsea plots
gsea:
	make -f gsea.mk IN=${IN} OUT=${OUT} ONT=${ONT} PREFIX=${PREFIX} plots

#
# Start from functional analysis results
#

# Create all ora plots starting from gprofiler GO results
ora_plots:
	make -f ora.mk ORA=${ORA_IN} OUT=${OUT} ONT=${ONT} DIRECTION=none PREFIX=${PREFIX} plots

