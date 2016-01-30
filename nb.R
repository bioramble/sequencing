############################################################################
# Bioramble
# Why sequencing data is modeled as negative binomial
# by Jesse Lipp
# created: Jan 30, 2016
############################################################################

# --------------------------------------------------------------------------
# Set up environment
# --------------------------------------------------------------------------
# clean-up
rm(list = ls())

# libraries
if (!require(edgeR)) {
  install.packages("edgeR")
  library(edgeR)
}

# --------------------------------------------------------------------------
# Mean-Variance plot
# --------------------------------------------------------------------------
# download bottomly count table from ReCount project
bottomly <- read.table("http://bowtie-bio.sourceforge.net/recount/countTables/bottomly_count_table.txt", row.names = 1, header = TRUE)
# download experimental information
experiment <- read.table("http://bowtie-bio.sourceforge.net/recount/phenotypeTables/bottomly_phenodata.txt", header = TRUE)
# counts of genes that are consistently expressed in all samples
expressed <- apply(bottomly, 1, function(r) all(r > 0))
counts <- bottomly[expressed, ]

# plot mean variance plot using function provided by edgeR
d <- DGEList(counts, group = experiment$strain)
d <- estimateCommonDisp(d)
plotMeanVar(d, show.raw.vars = TRUE, show.ave.raw.vars = FALSE, NBline = TRUE)
legend("topleft", legend = c("Poisson", "Negative Binomial"), fill = c("black", "dodgerblue"), border = c("black", "dodgerblue"), bty = "n")
