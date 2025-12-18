# Differential Gene Expression Analysis: Lung Cancer (Tumor vs Normal)

This repository contains an end-to-end differential gene expression (DGE) analysis comparing tumor vs. normal lung tissue samples. The analysis identifies significantly altered genes, explores biological pathways, and constructs gene–gene interaction networks to highlight potential mechanisms involved in cancer progression.

---

## Objective

To identify differentially expressed genes (DEGs) between tumor and normal samples and perform downstream analyses, including:

- Statistical modeling  
- Visualization  
- Functional enrichment (GO / KEGG)  
- Network analysis (STRING + Cytoscape)  

The goal is to understand which genes change significantly and what biological pathways they may influence.

---

## Project Overview

This project uses RNA-seq gene expression data and applies:

- `limma` for differential expression  
- `ggplot2` for visualization  
- `clusterProfiler` for GO / KEGG enrichment  
- STRING and Cytoscape for network analysis  

The workflow progresses from raw expression matrices → DEGs → pathway-level insights → network interactions.

---

## Repository Structure

    LungCancer_RNASeq_DEAnalysis/
    │
    ├── 01_data/
    │   ├── DEGs.csv
    │   ├── significant_genes.csv
    │   └── DEG_list_for_STRING.csv
    │
    ├── 02_scripts/
    │   ├── 01_data_loading.R
    │   ├── 02_data_preprocessing.R
    │   ├── 03_differential_expression_analysis.R
    │   ├── 04_gene_enrichment_analysis.R
    │   └── 05_network_analysis.R
    │
    ├── 03_results/
    │   └── visualizations.pdf
    │
    └── README.md

This structure reflects the chronological steps of the computational workflow.

---

## Steps in the Analysis

### 1. Data Loading & Inspection

**Script:** `01_data_loading.R`

- Imports the RNA-seq dataset (GEO)  
- Performs initial quality control  
- Checks sample annotations  

---

### 2. Data Preprocessing

**Script:** `02_data_preprocessing.R`

- Filters low-expression genes  
- Normalizes expression values  
- Maps gene identifiers for annotation  

Preprocessing ensures statistical reliability for limma modeling.

---

### 3. Differential Expression Analysis (limma)

**Script:** `03_differential_expression_analysis.R`

- Builds design matrix (Tumor vs Normal)  
- Fits linear model using `limma`  
- Computes log2 fold changes  
- Adjusts p-values using the Benjamini–Hochberg method  

Exports:

- `DEGs.csv`  
- `significant_genes.csv` (adjusted p-value < 0.05)  

Visualizations included:

- Volcano plot  
- Bar plot of top genes  
- Scatter plot of expression trends  

(All stored in `03_results/visualizations.pdf`)

---

### 4. Functional Enrichment Analysis

**Script:** `04_gene_enrichment_analysis.R`

Using `clusterProfiler`, the analysis includes:

- GO Biological Process enrichment  
- KEGG pathway enrichment  
- Visualization using dotplots and barplots  

Key enriched pathways include:

- Cytoskeleton organization  
- Immune-related processes  
- Glial differentiation  

These highlight potential mechanisms influencing lung cancer biology.

---

### 5. Network Analysis (STRING + Cytoscape)

**Script:** `05_network_analysis.R`

- Prepares DEG list with log2 fold-change values  
- Creates `DEG_list_for_STRING.csv`  
- Used in STRING to retrieve protein–protein interactions  
- Can be imported into Cytoscape for network visualization  

Network analysis helps reveal central hub genes and interaction modules.

---

## Data Description

The dataset contains RNA-seq gene expression values for tumor and normal lung tissue samples.

For each gene:

- Log2 fold change  
- P-value  
- Adjusted p-value (FDR)  

Significant DEGs (adjusted p-value < 0.05) are used for enrichment and network analysis.

---

## Key Files

- `DEGs.csv` — Full list of DEGs (tumor vs normal)  
- `significant_genes.csv` — DEGs with adjusted p-value < 0.05  
- `DEG_list_for_STRING.csv` — Input file for STRING / Cytoscape  
- `visualizations.pdf` — Volcano, bar plot, and scatter plot  

---

## Tools & Packages Used

### R Packages
- limma  
- ggplot2  
- clusterProfiler  
- org.Hs.eg.db  
- enrichplot  

### Network Tools
- STRING database  
- Cytoscape  

---

## Conclusion

This analysis identifies key genes that differ between tumor and normal lung samples and links them to biological pathways relevant to cancer progression. Enrichment results suggest involvement in:

- Cytoskeleton regulation  
- Immune system activation  
- Differentiation-related pathways  

These findings may serve as a starting point for deeper biological investigation.

---

## References

- limma documentation  
- clusterProfiler manual  
- Cytoscape documentation
