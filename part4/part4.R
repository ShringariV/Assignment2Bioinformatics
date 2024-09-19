# Set the CRAN mirror
options(repos = c(CRAN = "https://cran.rstudio.com/"))

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
# Use BiocManager to install ComplexHeatmap if it’s not installed
if (!("ComplexHeatmap" %in% installed.packages())) {
  BiocManager::install("ComplexHeatmap")
}
# Install RColorBrewer from CRAN if it’s not installed
if (!("RColorBrewer" %in% installed.packages())) {
  install.packages("RColorBrewer")
}
# Install the 'readr' package if it's not already installed
if (!("readr" %in% installed.packages())) {
  install.packages("readr")
}
# Install the 'dplyr' package if it's not already installed
if (!("dplyr" %in% installed.packages())) {
  install.packages("dplyr")
}
# Install the 'magrittr' package if it's not already installed
if (!("magrittr" %in% installed.packages())) {
  install.packages("magrittr")
}
# Install the 'tibble' package if it's not already installed
if (!("tibble" %in% installed.packages())) {
  install.packages("tibble")
}

# Load the required libraries
library(ComplexHeatmap)
library(RColorBrewer)
library(readr)
library(dplyr)
library(magrittr) # For %>%
library(tibble) # For rownames_to_column and column_to_rownames

# Path to Part 3 folder
part3_folder <- "../part3"
# Path to data folder
data_folder <- "../data"
setwd("/Users/joelalomafernandez/Desktop/Assignment2Bio/Assignment2Bioinformatics/part4")

# Load the DESeq2 results from the Part 3 folder
deseq_df <- readr::read_tsv(file.path(part3_folder, "ERP107715_diff_expr_results.tsv"), show_col_types = FALSE)

# Filter the significant genes from the DESeq2 results (padj < 0.05)
significant_genes <- deseq_df %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(padj)

# Take the top 50 most significant genes for the heatmap
top_genes <- significant_genes %>%
  dplyr::slice(1:50) %>%
  dplyr::pull(Gene)

# Load the expression matrix
expression_df <- readr::read_tsv(file.path(data_folder, "ERP107715.tsv"), show_col_types = FALSE)

# Extract the gene expression matrix for the top 50 genes
top_gene_expression <- expression_df %>%
  dplyr::filter(Gene %in% top_genes) %>%
  tibble::column_to_rownames(var = "Gene")  # Set Gene as row names

# Normalize the expression data (optional but improves visualization)
# Here, we're using z-score normalization for each gene
top_gene_expression_scaled <- t(apply(top_gene_expression, 1, scale))

# Assign proper row names and column names to the scaled data
rownames(top_gene_expression_scaled) <- top_genes
colnames(top_gene_expression_scaled) <- colnames(top_gene_expression)

# Define a color palette for the heatmap
col_palette <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)

# Load the metadata file
metadata_file <- file.path(data_folder, "metadata_ERP107715.tsv") # Adjust file path as needed
metadata_df <- readr::read_tsv(metadata_file, show_col_types = FALSE)

# Extract relevant columns for annotations
# Assuming 'refinebio_accession_code' matches the column names in your expression data
annotation_df <- metadata_df %>%
  dplyr::select(refinebio_accession_code, refinebio_treatment) %>%
  dplyr::rename(Sample = refinebio_accession_code, Treatment = refinebio_treatment)

# Create a named vector for annotations
annotation_vector <- setNames(annotation_df$Sex, annotation_df$Sample)

# Create the annotation object for the heatmap
annotation <- ComplexHeatmap::HeatmapAnnotation(
  Sex = annotation_vector,
  col = list(Sex = c("female" = "pink", "male" = "blue")) # Adjust colors as needed
)

# Generate the heatmap
heatmap_plot <- ComplexHeatmap::Heatmap(
  top_gene_expression_scaled, 
  name = "Expression", # Title for the heatmap legend
  col = col_palette,   # Color palette for the heatmap
  show_row_names = TRUE, # Show the gene names (rows)
  show_column_names = TRUE, # Show sample names (columns)
  row_names_gp = gpar(fontsize = 10), # Control gene name font size
  column_names_gp = gpar(fontsize = 10), # Control column name font size
  clustering_method_rows = "complete", # Hierarchical clustering method for rows
  clustering_method_columns = "complete", # Hierarchical clustering method for columns
  heatmap_legend_param = list(
    title = "Value", # Title of the legend
    title_gp = gpar(fontsize = 10), # Font size of the title
    labels_gp = gpar(fontsize = 10), # Font size of the labels
    legend_height = unit(5, "cm") # Height of the legend
  ),
  top_annotation = annotation # Add the sidebar annotation
)

# Define the file path to save the heatmap
file_path <- file.path(getwd(), "ERP107715_Heatmap.png")

# Save the heatmap as a PNG image
png(filename = file_path, width = 1200, height = 1000) # Adjust the width and height as needed
draw(heatmap_plot)
dev.off()