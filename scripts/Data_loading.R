# Set working directory
# Load required libraries
library(limma)
library(tidyverse)
# Load the data (expression matrix)
series_matrix <- read.table("GSE19804_series_matrix.txt.gz", 
                            header = TRUE, sep = "\t", 
                            comment.char = "!", 
                            stringsAsFactors = FALSE)
# Preview the first few rows of the data
head(series_matrix)

# Read the metadata (header lines)
meta_lines <- readLines("GSE19804_series_matrix.txt.gz")
# Search for sample characteristics (Normal vs Tumor info)
grep("!Sample_characteristics_ch1", meta_lines, value = TRUE)
