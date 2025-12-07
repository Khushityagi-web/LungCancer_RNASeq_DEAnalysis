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
library(ggplot2)
# Volcano plot code
ggplot(top_genes, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("red", "gray")) +
  theme_minimal() +
  ggtitle("Volcano Plot: Tumor vs Normal") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted p-value")
ggsave("figures/volcano_plot.png")
ggplot(top_genes_filtered, aes(x = reorder(Gene, logFC), y = logFC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  xlab("Genes") +
  ylab("Log2 Fold Change")
ggsave("figures/bar_plot.png")
plot(top_genes$logFC, -log10(top_genes$adj.P.Val), pch = 20, 
     xlab = "Log2 Fold Change", ylab = "-Log10 Adjusted P-value")
