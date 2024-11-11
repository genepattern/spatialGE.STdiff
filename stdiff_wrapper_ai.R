#!/usr/bin/env Rscript

# Load necessary libraries
library(optparse)
library(spatialGE) # Assuming STdiff is available as a package
library(ggplot2) # For plotting

# Define command line options
option_list <- list(
  make_option(c("-x", "--x"), type="character", default=NULL,
              help="STlist object containing spatial transcriptomics data", metavar="character"),
  make_option(c("-s", "--samples"), type="character", default=NULL,
              help="Comma-separated list of sample indices or names to include in DE tests", metavar="character"),
  make_option(c("-a", "--annot"), type="character", default=NULL,
              help="Column name in x@spatial_meta containing groups/clusters to be tested", metavar="character"),
  make_option(c("-w", "--w"), type="numeric", default=NULL,
              help="Spatial weight used in STclust. Required if annot is empty", metavar="numeric"),
  make_option(c("-k", "--k"), type="character", default=NULL,
              help="k value used in STclust, or 'dtc' for dynamicTreeCut clusters. Required if annot is empty", metavar="character"),
  make_option(c("-d", "--deepSplit"), type="logical", default=NULL,
              help="deepSplit value if used in STclust. Required if k='dtc'", metavar="logical"),
  make_option(c("-t", "--topgenes"), type="integer", default=5000,
              help="Number of top variable genes to select from each sample based on variance", metavar="integer"),
  make_option(c("-p", "--pval_thr"), type="numeric", default=0.05,
              help="Cut-off of adjusted p-values to define differentially expressed genes", metavar="numeric"),
  make_option(c("-j", "--pval_adj"), type="character", default="fdr",
              help="Method to adjust p-values. Defaults to 'FDR'", metavar="character"),
  make_option(c("-g", "--sp_topgenes"), type="numeric", default=0.2,
              help="Proportion of DE genes from non-spatial models to use in spatial models", metavar="numeric"),
  make_option(c("-c", "--clusters"), type="character", default=NULL,
              help="Comma-separated list of cluster names to test DE genes", metavar="character"),
  make_option(c("-r", "--pairwise"), type="logical", default=FALSE,
              help="Whether to carry tests on a pairwise manner", metavar="logical"),
  make_option(c("-v", "--verbose"), type="integer", default=1,
              help="Output progress indicators. 0 for no text, 1 for basic, 2 for detailed", metavar="integer"),
  make_option(c("-n", "--cores"), type="integer", default=NULL,
              help="Number of cores to use in parallelization", metavar="integer"),
  make_option(c("-P", "--plot"), type="logical", default=FALSE,
              help="Option to plot intermediate results", metavar="logical")
)

# Parse command line options
opt <- parse_args(OptionParser(option_list=option_list))

# Convert samples and clusters from comma-separated strings to vectors
if (!is.null(opt$samples)) {
  opt$samples <- unlist(strsplit(opt$samples, ","))
}
if (!is.null(opt$clusters)) {
  opt$clusters <- unlist(strsplit(opt$clusters, ","))
}

# Load the STlist object
# Assuming the STlist object is saved in an RData file
data <- readRDS(opt$x)

# Run the STdiff function
result <- STdiff(
  x = data,
  samples = opt$samples,
  annot = opt$annot,
  w = opt$w,
  k = opt$k,
  deepSplit = opt$deepSplit,
  topgenes = opt$topgenes,
  pval_thr = opt$pval_thr,
  pval_adj = opt$pval_adj,
  sp_topgenes = opt$sp_topgenes,
  clusters = opt$clusters,
  pairwise = opt$pairwise,
  verbose = opt$verbose,
  cores = opt$cores
)

# Plot intermediate results if requested
if (opt$plot) {
  # Example plot: Number of DE genes per sample
  de_counts <- sapply(result, nrow)
  ggplot(data.frame(Sample=names(de_counts), DE_Genes=de_counts), aes(x=Sample, y=DE_Genes)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    labs(title="Number of Differentially Expressed Genes per Sample", x="Sample", y="DE Genes")
}

save(result, file="/STdiff_results.RData")

### R Command Line to Run the Wrapper Script
# Rscript STdiff_wrapper.R --x "path/to/STlist.RData" --samples "1,2,3" --annot "cluster_column" --topgenes 5000 --pval_thr 0.05 --pval_adj "fdr" --sp_topgenes 0.2 --verbose 1 --plot TRUE

### Explanation
# - The script uses `optparse` to define and parse command-line options.
# - It loads the `STlist` object from a specified file.
# - It runs the `STdiff` function with the provided parameters.
# - If the `--plot` option is set to `TRUE`, it generates a plot of the number of differentially expressed genes per sample.
# - The results are saved to an RData file for further analysis.