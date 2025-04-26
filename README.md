## Differential Gene Expression (DGE) Analysis: Tumor vs. Normal
# Objective
The objective of this project is to identify the differentially expressed genes (DEGs) between tumor and normal samples and conduct downstream analysis, including visualization, functional enrichment, and network analysis. The project aims to understand which genes show significant expression changes between the two groups and their potential biological significance in cancer-related pathways.
# Project Overview
This project uses RNA-Seq data to perform differential expression analysis and subsequent functional enrichment of the DEGs. We employ tools and packages such as limma, ggplot2, clusterProfiler, and Cytoscape to perform the analysis and visualize the results.
# Steps in the Analysis:
Steps in the Analysis:
Data Preprocessing: Filter out low-expression genes and normalize the data.

Differential Expression Analysis: Using limma to identify DEGs between tumor and normal samples.

Visualization:

Volcano plot

Bar plot

Scatter plot

Functional Enrichment Analysis: Perform GO and KEGG enrichment using clusterProfiler.

Network Analysis: Identify key gene-gene interactions using Cytoscape and STRING database.
# Data Description
The data used in this analysis consists of RNA-Seq expression levels of genes across tumor and normal samples. The data was processed to obtain log2 fold changes and adjusted p-values for each gene. The significant DEGs (adjusted p-value < 0.05) were further analyzed for functional enrichment and gene-gene interactions.
# Code Files
data_loading.R: Script for importing and inspecting the GEO dataset.

data_preprocessing.R: Script for preprocessing microarray data, including normalization and annotation.

differential_expression_analysis.R: Performs differential expression analysis using limma.

gene_enrichment_analysis.R: Conducts GO/KEGG enrichment using clusterProfiler.

network_analysis.R: Prepares input for STRING and Cytoscape to analyze DEG networks.
# Key Files
DEGs.csv: A list of differentially expressed genes (DEGs) between tumor and normal samples.

significant_genes.csv: Exported list of significant DEGs (adjusted p-value < 0.05).

DEG_list_for_STRING.csv: List of DEGs with log2 fold change values for network analysis.

functional_enrichment_results.csv: Results from Gene Ontology (GO) enrichment analysis.

visualizations: Folder containing the generated plots (volcano plot, bar plot, scatter plot, GO enrichment plots).
Visualizations
The following visualizations are included:

Volcano Plot

Bar Plot

Scatter Plot

GO Enrichment Bar and Dot Plot
# Conclusion
This analysis identifies key genes that are differentially expressed between tumor and normal samples and provides insights into the biological processes these genes are involved in. Functional enrichment analysis points to significant pathways related to cytoskeleton regulation, immune response, and glial differentiation, which may have implications in cancer research.
# References
Limma Package Documentation: https://bioconductor.org/packages/release/bioc/html/limma.html

ClusterProfiler Documentation: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html

Cytoscape: https://cytoscape.org/
