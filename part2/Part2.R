
# Define the file path to the data directory
# Replace with the path of the folder the files will be in
data_dir <- file.path("dataSRP", "")

# Declare the file path to the gene expression matrix file
# inside directory saved as `data_dir`
# Replace with the path to your dataset file
data_file <- file.path(data_dir, "SRP111402.tsv")

# Declare the file path to the metadata file
# inside the directory saved as `data_dir`
# Replace with the path to your metadata file
metadata_file <- file.path(data_dir, "metadata_SRP111402.tsv")

file.exists(data_file)
file.exists(metadata_file)

if (!("DESeq2" %in% installed.packages())) {
  # Install DESeq2
  BiocManager::install("DESeq2", update = FALSE)
}
# Attach the `DESeq2` library
library(DESeq2)

# Attach the `ggplot2` library for plotting
library(ggplot2)

# We will need this so we can use the pipe: %>%
library(magrittr)

# Set the seed so our results are reproducible:
set.seed(12345)
# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)
# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  # Tuck away the gene ID column as row names, leaving only numeric values
  tibble::column_to_rownames("Gene")

# Make the sure the columns (samples) are in the same order as the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

# convert the columns we will be using for annotation into factors
metadata <- metadata %>%
  dplyr::mutate(
    refinebio_specimen_part = factor(
      refinebio_specimen_part,
      levels = c("metastatic", "non-malignant")
    ),
    refinebio_title = as.factor(refinebio_title)
  )
# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)

# The `DESeqDataSetFromMatrix()` function needs the values to be integers
filtered_expression_df <- round(filtered_expression_df)
# Create a `DESeqDataSet` object
dds <- DESeqDataSetFromMatrix(
  countData = filtered_expression_df, # the counts values for all samples in our dataset
  colData = metadata, # annotation data for the samples in the counts data frame
  design = ~1 # Here we are not specifying a model
  # Replace with an appropriate design variable for your analysis
)

# Normalize and transform the data in the `DESeqDataSet` object
# using the `vst()` function from the `DESeq2` R package
dds_norm <- vst(dds)

pca_results <-
  plotPCA(
    dds_norm,
    intgroup = c("refinebio_specimen_part", "refinebio_title"),
    returnData = TRUE
  )

annotated_PCA_plot <- ggplot(
  pca_results,
  aes(
    x= PC1,
    y= PC2,
    color = refinebio_specimen_part,
    shape = refinebio_title
  )
) +
  geom_point()

annotated_PCA_plot
plots_dir <- "plots"
ggsave(
  file.path(plots_dir, "SRP111402_pca_plot.png"),
  # Replace with a file name relevant your plotted data
  plot = annotated_PCA_plot # the plot object that we want saved to file
)
