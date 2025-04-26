write.csv(sig_genes_mapped[, c("GeneSymbol", "logFC")], file = "data/DEG_list_for_STRING.csv", row.names = FALSE)
