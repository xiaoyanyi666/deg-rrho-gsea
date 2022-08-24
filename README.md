# deg-rrho-gsea

Description

Script.R: a workflow pipeline for differential gene expression analysis, rank-rank hypergeometric overlap analysis, and gene set enrichment analysis

Methods

1) The input to the script is tab-separated metadata.tsv 
2) The salmon quantification for each sample is quired
3) A homemade R package “rucdr” should be installed to run this script
4) Other required packages and databases referred to the script.R

Installation of “rucdr”

1) Download the package from the official depository https://github.com/antpiron/rucdr
2) sudo R CMD INSTALL rucdr_0.1.tar.gz for Debian/Ubuntu
3) install.packages(“rucdr_0.1.tar.gz”, repos= NULL, type=”source”) for R
3) library(“rucdr”)
