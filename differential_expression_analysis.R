# Load limma if not already loaded
library(limma)

# Prepare the expression data matrix
exprs_matrix <- as.matrix(expr_data[, -1])  # remove the ID column if it's not needed
rownames(exprs_matrix) <- expr_data[, 1]    # set gene IDs as row names

# Create design matrix
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "Tumor_vs_Normal")

# Fit the linear model
fit <- lmFit(exprs_matrix, design)

# Apply empirical Bayes smoothing
fit <- eBayes(fit)

# View top differentially expressed genes
top_genes <- topTable(fit, coef="Tumor_vs_Normal", adjust="fdr", number=20)

# Save the top genes to a CSV
write.csv(top_genes, "results/top_differentially_expressed_genes.csv")

# View results
head(top_genes)
