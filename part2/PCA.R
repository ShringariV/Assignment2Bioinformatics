results_dir <- file.path("results", "")

# Ensure the directory exists
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}
# Define the file path to the data directory
data_dir <- file.path("data", "")

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}


# Declare the file path to the gene expression matrix file
# inside directory saved as `data_dir`
# Replace with the path to your dataset file
data_file <- file.path(data_dir, "ERP107715.tsv")

# Declare the file path to the metadata file
# inside the directory saved as `data_dir`
# Replace with the path to your metadata file
metadata_file <- file.path(data_dir, "metadata_ERP107715.tsv")

#If exists
file.exists(data_file)
file.exists(metadata_file)

if(!("DESeq2" %in% installed.packages())){
  #Install DESeq2
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
    refinebio_sex = factor(
      refinebio_sex,
      # specify the possible levels in the order we want them to appear
      levels = c("female", "male")
    ),
    refinebio_disease = as.factor(refinebio_disease)
  )

# Define a minimum counts cutoff and filter the data to include
# only rows (genes) that have total counts above the cutoff
filtered_expression_df <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)

filtered_expression_df <- round(filtered_expression_df)
# Check if any values in the countData are negative
sum(filtered_expression_df < 0)
# Replace negative values with 0 (if raw counts are expected)
filtered_expression_df[filtered_expression_df < 0] <- 0
# Ensure count data is an integer matrix
filtered_expression_df <- round(filtered_expression_df)

table(metadata$refinebio_disease)

dds <- DESeqDataSetFromMatrix(
  countData = filtered_expression_df,
  colData = metadata,
  design = ~ 1
)

dds_norm <- vst(dds)
pca_results <-
  plotPCA(
    dds_norm,
    intgroup = c("refinebio_sex", "refinebio_disease"),
    returnData = TRUE # This argument tells R to return the PCA values
  )
# Plot using `ggplot()` function and save to an object
annotated_pca_plot <- ggplot(
  pca_results,
  aes(
    x = PC1,
    y = PC2,
    # plot points with different colors for each `refinebio_treatment` group
    color = refinebio_sex,
    # plot points with different shapes for each `refinebio_disease` group
    shape = refinebio_disease
  )
) +
  # Make a scatter plot
  geom_point()

# display annotated plot
annotated_pca_plot
# Save plot using `ggsave()` function
ggsave(
  file.path(plots_dir, "ERP107715_pca_plot.png"),
  # Replace with a file name relevant your plotted data
  plot = annotated_pca_plot # the plot object that we want saved to file
)
