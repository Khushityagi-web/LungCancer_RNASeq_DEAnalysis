write.csv(sig_genes, "data/significant_genes.csv", row.names = TRUE)
# GO/KEGG enrichment analysis code
barplot(go_bp_all, showCategory = 10, title = "Top GO BP Terms (Unfiltered)")
ggsave("figures/enrichment_dotplot.png")
library(AnnotationDbi)
library(hgu133plus2.db)
gene_symbols <- mapIds(hgu133plus2.db, keys = clean_probe_ids, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
sig_genes$GeneSymbol <- gene_symbols
write.csv(sig_genes, "data/significant_genes.csv", row.names = TRUE)


