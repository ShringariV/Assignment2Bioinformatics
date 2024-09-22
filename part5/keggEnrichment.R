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
if (!("biomaRt" %in% installed.packages())) {
  BiocManager::install("biomaRt")
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(dplyr)
library(ggplot2)

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
deg_genes_symbols <- unique(na.omit(deseq_with_symbols$external_gene_name))

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(deg_genes_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Get unique Entrez IDs
deg_genes_entrez <- unique(entrez_ids$ENTREZID)

# Perform KEGG pathway analysis
kegg_results <- enrichKEGG(
  gene = deg_genes_entrez,
  organism = "hsa",  # Homo sapiens
  pAdjustMethod = "BH",  # Benjamini-Hochberg adjustment
  qvalueCutoff = 0.05  # Q-value cutoff
)

# View the top results
print(head(kegg_results))

# Visualize the results
# Bar plot of the top 10 KEGG pathways
barplot(kegg_results, showCategory = 10) +
  ggtitle("Top 10 KEGG Enriched Pathways")

# Dot plot of KEGG enrichment results
dotplot(kegg_results, showCategory = 10) +
  ggtitle("KEGG Enrichment Dot Plot")

# Save results to a CSV file
write.csv(as.data.frame(kegg_results), file = "kegg_enrichment_results.csv", row.names = FALSE)