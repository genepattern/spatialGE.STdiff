# Description
This script performs differential expression analysis on spatial transcriptomics data using the `spatialGE` package. It leverages the `STdiff` function from `spatialGE` to identify differentially expressed genes based on spatial and non-spatial models. The output is an `.RDS` file containing results for further downstream analysis and visualization.

# Module Details
- **Authors:** Edwin Huang, Aditya Kakarla, GenePattern LLM Pipeline  
- **Categories:** Spatial transcriptomics
- **Source Repo:** [GitHub: spatialGE](https://github.com/FridleyLab/spatialGE)  
- **Contact:** [GenePattern Help](https://groups.google.com/u/1/g/genepattern-help)  
- **Programming Language:** R  

# Input Files
| Info              | Type     |
|------------------|---------|
| input.RDS       | Required |

# Output Files
| Info                          | Type     |
|-----------------------------|---------|
| STdiff_results.rds          | Required |

# Parameters
| Parameter | Name | Description | Default Value | Type |
|-----------|------|-------------|---------------|------|
| --x | -x | STlist object containing spatial transcriptomics data | None | str |
| --s | --samples | Comma-separated list of sample indices or names to include in DE tests | None | str |
| --a | --annot | Column name in `x@spatial_meta` containing groups/clusters to be tested | None | str |
| --w | --w | Spatial weight used in STclust. Required if `annot` is empty | None | float |
| --k | --k | k value used in STclust, or 'dtc' for dynamicTreeCut clusters. Required if `annot` is empty | None | str |
| --d | --deepSplit | deepSplit value if used in STclust. Required if `k='dtc'` | None | bool |
| --t | --topgenes | Number of top variable genes to select from each sample based on variance | 5000 | int |
| --p | --pval_thr | Cut-off of adjusted p-values to define differentially expressed genes | 0.05 | float |
| --j | --pval_adj | Method to adjust p-values (e.g., 'fdr') | 'fdr' | str |
| --g | --sp_topgenes | Proportion of DE genes from non-spatial models to use in spatial models | 0.2 | float |
| --c | --clusters | Comma-separated list of cluster names to test DE genes | None | str |
| --r | --pairwise | Whether to carry tests on a pairwise manner | FALSE | bool |
| --v | --verbose | Output progress indicators: 0 (no text), 1 (basic), 2 (detailed) | 1 | int |
| --n | --cores | Number of cores to use in parallelization | None | int |
| --P | --plot | Option to plot intermediate results | FALSE | bool |