# Load Required Libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!("clusterProfiler" %in% installed.packages())) {
  BiocManager::install("clusterProfiler")
}
if (!("org.Hs.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Hs.eg.db")
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(biomaRt)

result_dir <- file.path("results","")
# Define the file path to the data directory
data_dir <- file.path("data", "")

# Declare the file path to the gene expression matrix file
data_file <- file.path(data_dir, "ERP107715.tsv")
result_file <- file.path(result_dir, "ERP107715_diff_expr_results.tsv")

deseq_df <- readr::read_tsv(result_file, show_col_types = FALSE)

# Filter for significant genes (padj < 0.05)
significant_genes <- deseq_df %>%
  filter(padj < 0.05) %>%
  pull(Gene)  # Assuming "Gene" column contains Ensembl IDs

# Remove NA values
significant_genes <- na.omit(significant_genes)

# Map Ensembl IDs to Gene Symbols
# Use biomaRt to get gene symbols
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene symbols
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = significant_genes,
  mart = ensembl
)

# Merge with DESeq2 results to get symbols for significant genes
deseq_with_symbols <- deseq_df %>%
  left_join(gene_map, by = c("Gene" = "ensembl_gene_id"))

# Get unique gene symbols for enrichment analysis
deg_genes <- unique(na.omit(deseq_with_symbols$external_gene_name))

# Perform GO enrichment analysis
go_results <- enrichGO(
  gene = deg_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",  # Use "SYMBOL" for gene symbols
  ont = "BP",  # Biological process
  pAdjustMethod = "BH",  # Benjamini-Hochberg adjustment
  qvalueCutoff = 0.05,  # Q-value cutoff
  readable = TRUE  # Convert gene IDs to readable names
)

# View the top results
print(head(go_results))

# Visualize the results
# Bar plot of the top 10 GO terms
barplot(go_results, showCategory = 10) +
  ggtitle("Top 10 GO Enriched Terms")

# Dot plot of GO enrichment results
dotplot(go_results, showCategory = 10) +
  ggtitle("GO Enrichment Dot Plot")

# Save results to a CSV file
write.csv(as.data.frame(go_results), file = "go_enrichment_results.csv", row.names = FALSE)