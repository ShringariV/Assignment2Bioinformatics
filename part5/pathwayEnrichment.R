library(fgsea)
library(org.Hs.eg.db)
library(msigdbr)
library(dplyr)
library(tibble)
library(biomaRt)

result_dir <- file.path("results","")
# Define the file path to the data directory
data_dir <- file.path("data", "")

# Declare the file path to the gene expression matrix file
data_file <- file.path(data_dir, "ERP107715.tsv")
result_file <- file.path(result_dir, "ERP107715_diff_expr_results.tsv")
# Load the DESeq2 results from the Part 3 folder
deseq_df <- readr::read_tsv(result_file, show_col_types = FALSE)

ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

ensemblToGene <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = deseq_df$Gene,
  mart = ensembl
)

# Merge the gene symbol information with your DESeq2 results
deseq_with_symbols <- deseq_df %>%
  dplyr::left_join(ensemblToGene, by = c("Gene" = "ensembl_gene_id")) %>%
  dplyr::filter(hgnc_symbol != "")  # Remove rows without mapped symbols

# Step 2: Prepare the ranked gene list for GSEA
gene_ranking <- deseq_with_symbols %>%
  dplyr::filter(!is.na(log2FoldChange)) %>%  # Remove rows with NA in log2foldchange
  dplyr::arrange(desc(log2FoldChange)) %>%   # Sort genes by log2 fold change
  dplyr::mutate(gene_rank = log2FoldChange) %>%  # Create a ranking column
  dplyr::select(hgnc_symbol, gene_rank) %>%  # Select gene symbol and log2foldchange (gene_rank)
  tibble::deframe()  # Convert to named vector

# Step 3: Load MSigDB gene sets for human
msigdb_gene_sets <- msigdbr(species = "human", category = "H")

# Convert msigdbr result to a list of pathways for fgsea
pathways <- split(msigdb_gene_sets$gene_symbol, msigdb_gene_sets$gs_name)

# Step 4: Perform GSEA using fgsea
fgsea_results <- fgsea(
  pathways = pathways,
  stats = gene_ranking,
  minSize = 10,
  maxSize = 500,
  nperm = 1000  # Number of permutations
)

# Sort by adjusted p-value
fgsea_results <- fgsea_results %>% arrange(padj)

# Print the top GSEA results
print(head(fgsea_results))

# Step 5: Plot the top pathway (if any significant)
if (nrow(fgsea_results) > 0) {
  topPathway <- fgsea_results$pathway[1]  # Select the top enriched pathway
  plotEnrichment(pathways[[topPathway]], gene_ranking) +
    ggplot2::ggtitle(paste("Top Pathway:", topPathway))
} else {
  cat("No significantly enriched pathways found.\n")
}
