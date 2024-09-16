# Define the directory where the results will be saved
results_dir <- file.path("results", "")

# Ensure the directory exists
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}
# Define the file path to the data directory
data_dir <- file.path("data", "")

# Declare the file path to the gene expression matrix file
data_file <- file.path(data_dir, "ERP1107715.tsv")

# Declare the file path to the metadata file
metadata_file <- file.path(data_dir, "metadata_ERP1107715.tsv")

# Check if files exist
file.exists(data_file)
file.exists(metadata_file)

# Install org.Hs.eg.db if not already installed
if (!("org.Hs.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

library(org.Hs.eg.db)
library(magrittr)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)

# Read in metadata TSV file
metadata <- read_tsv(metadata_file)

# Read in data TSV file
expression_df <- read_tsv(data_file) %>%
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

# Bring back the "Gene" column in preparation for mapping
expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

# Map Ensembl IDs to their associated Gene Names
mapped_list <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = expression_df$Gene,
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "GENENAME", # Map to Gene Names instead of ENTREZID
  multiVals = "list"
)

# Preview the mapped list
head(mapped_list)

# Convert the mapped list to a data frame
mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "GeneName") %>%
  tidyr::unnest(cols = GeneName)

# Show the distribution of Gene Names
summary(as.factor(mapped_df$GeneName), maxsum = 10)

# Count the number of times each Ensembl ID appears
multi_mapped <- mapped_df %>%
  dplyr::count(Ensembl, name = "gene_name_count") %>%
  dplyr::arrange(desc(gene_name_count))

# Preview the first 6 rows of multi_mapped
head(multi_mapped)

# Collapse Gene Names into one column
collapsed_mapped_df <- mapped_df %>%
  dplyr::group_by(Ensembl) %>%
  dplyr::summarize(all_gene_names = paste(GeneName, collapse = ";"))

# Preview rows with multiple Gene Names
collapsed_mapped_df %>%
  dplyr::filter(str_detect(all_gene_names, ";")) %>%
  head()

# Create final mapped data frame with the first Gene Name and all Gene Names
final_mapped_df <- data.frame(
  "first_mapped_gene_name" = mapIds(
    org.Hs.eg.db, # Replace with annotation package for your organism
    keys = expression_df$Gene,
    keytype = "ENSEMBL", # Replace with the gene identifiers used in your data
    column = "GENENAME", # Map to Gene Names
    multiVals = "first" # Keep only the first mapped value for each Ensembl ID
  )
) %>%
  tibble::rownames_to_column("Ensembl") %>%
  dplyr::inner_join(collapsed_mapped_df, by = "Ensembl") %>%
  dplyr::inner_join(expression_df, by = c("Ensembl" = "Gene"))

# Preview the final mapped data frame with multiple Gene Names
final_mapped_df %>%
  dplyr::filter(str_detect(all_gene_names, ";")) %>%
  head()
# Write mapped and annotated data frame to output file
readr::write_tsv(final_mapped_df, file.path(
  results_dir,
  "GSE140684_GeneNames.tsv" # Replace with a relevant output file name
))
