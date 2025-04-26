# Extract sample labels
lung_cancer_samples <- grep("lung cancer", meta_lines)
normal_adjacent_samples <- grep("paired normal adjacent", meta_lines)

# Extract corresponding data
lung_cancer_data <- meta_lines[lung_cancer_samples]
normal_adjacent_data <- meta_lines[normal_adjacent_samples]

# Process the expression data
lines <- readLines("GSE19804_series_matrix.txt")
start <- grep("!series_matrix_table_begin", lines) + 1
end <- grep("!series_matrix_table_end", lines) - 1
expr_lines <- lines[start:end]
expr_data <- read.table(text = expr_lines, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)

# Organize metadata for analysis
metadata_lines <- grep("!Sample_characteristics_ch1", lines, value = TRUE)
sample_labels <- sapply(strsplit(metadata_lines, "\t"), function(x) x[-1])
sample_labels <- unlist(sample_labels)

tissue_labels <- sample_labels[seq(1, length(sample_labels), by = 4)]
tissue_labels_clean <- gsub("\"", "", tissue_labels)

# Group into Tumor vs Normal
group <- ifelse(tissue_labels_clean == "tissue: lung cancer", "Tumor", "Normal")
group <- factor(group)
table(group)
