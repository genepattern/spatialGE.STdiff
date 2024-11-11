# options(timeout=9999999) # To avoid R closing connection with GitHub
# devtools::install_github("fridleylab/spatialGE")
library(spatialGE)

# Set directory path for melanoma data
data_files <- './data/melanoma_thrane/'

# List files for counts, mapping, and clinical data
count_files <- list.files(data_files, full.names = TRUE, pattern = 'counts')
coord_files <- list.files(data_files, full.names = TRUE, pattern = 'mapping')
clin_file <- list.files(data_files, full.names = TRUE, pattern = 'clinical')

# Load and save melanoma data
melanoma <- STlist(rnacounts=count_files, spotcoords=coord_files, samples=clin_file)
melanoma <- transform_data(melanoma)
melanoma <- STclust(melanoma, 
                    ks='dtc', 
                    ws=0.025)
saveRDS(melanoma, file = "melanoma_data.rds", version = 2)