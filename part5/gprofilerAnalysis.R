library(gprofiler2)
library(dplyr)
library(biomaRt)
library(tibble)
library(ggplot2)

result_dir <- file.path("results","")
data_dir <- file.path("data", "")

# Declare the file path to the gene expression matrix file
data_file <- file.path(data_dir, "ERP107715.tsv")
result_file <- file.path(result_dir, "ERP107715_diff_expr_results.tsv")
# Load the DESeq2 results from the Part 3 folder
deseq_df <- readr::read_tsv(result_file, show_col_types = FALSE)

# Use biomaRt to map Ensembl IDs to gene symbols
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

# Prepare the ranked gene list for g:Profiler
gene_ranking <- deseq_with_symbols %>%
  dplyr::filter(!is.na(log2FoldChange)) %>%  # Remove rows with NA in log2foldchange
  dplyr::arrange(desc(log2FoldChange)) %>%   # Sort genes by log2 fold change
  dplyr::mutate(gene_rank = log2FoldChange) %>%  # Create a ranking column
  dplyr::select(hgnc_symbol, gene_rank) %>%  # Select gene symbol and log2foldchange (gene_rank)
  tibble::deframe()  # Convert to named vector

# Convert ranked genes into a vector of gene symbols
ranked_genes <- names(gene_ranking)

# Perform enrichment analysis using gprofiler2's gost function
gost_res <- gost(
  query = ranked_genes,          # Ranked gene list
  organism = "hsapiens",          # Human
  ordered_query = TRUE,           # Indicate that the query is ranked
  significant = TRUE,             # Only keep significant results
  sources = c("GO:BP", "KEGG")    # You can specify different sources like GO terms, KEGG pathways, etc.
)
# Define file paths for saving results and plots
csv_output_path <- file.path(result_dir, "enrichment results/gprofiler_results.csv")
plot_output_path <- file.path(result_dir, "enrichment results/plots/gprofiler_enrichment_plot.png")

library(tidyr)

# Check if there are significant results and print them
if (!is.null(gost_res$result)) {
  # Convert any list columns to strings so that they can be written to a CSV
  flat_results <- gost_res$result %>%
    dplyr::mutate(across(where(is.list), ~ sapply(., toString)))  # Convert list columns to comma-separated strings

  # Save the flattened results to a CSV file
  write.csv(flat_results, csv_output_path, row.names = FALSE)
  cat("Results saved to:", csv_output_path, "\n")

  # Print the top results
  print(head(flat_results))
} else {
  cat("No significant enrichment found.\n")
}

# Visualize the top pathways using a bar plot (if results exist)
if (!is.null(gost_res$result) && nrow(gost_res$result) > 0) {
  # Create the plot
  enrichment_plot <- gostplot(gost_res, capped = TRUE, interactive = FALSE) +
    ggplot2::ggtitle("Top Enriched Pathways")

  # Save the plot as an image (png)
  ggsave(plot_output_path, plot = enrichment_plot, width = 10, height = 8, dpi = 300)
  cat("Plot saved to:", plot_output_path, "\n")

  # Display the plot
  print(enrichment_plot)
} else {
  cat("No enriched pathways found to plot.\n")
}