# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Define the file path to the data directory
data_dir <- file.path("dataGSE", "")

# Declare the file path to the gene expression matrix file
data_file <- file.path(data_dir, "GSE140684.tsv")

# Declare the file path to the metadata file
metadata_file <- file.path(data_dir, "metadata_GSE140684.tsv")

# Check if the gene expression matrix file is at the path stored in `data_file`
file.exists(data_file)
# Check if the metadata file is at the file path stored in `metadata_file`
file.exists(metadata_file)

if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("apeglm", update = FALSE)
}

# Attach the DESeq2 library
library(DESeq2)

# Attach the ggplot2 library for plotting
library(ggplot2)

# We will need this so we can use the pipe: %>%
library(magrittr)

set.seed(12345)

# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)

# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

head(metadata$refinebio_title)

metadata <- metadata %>%
  # Let's get the Stelara and placebo status from this variable
  dplyr::mutate(treatment_status = dplyr::case_when(
    stringr::str_detect(refinebio_title, "Stelara") ~ "Stelara",
    stringr::str_detect(refinebio_title, "Placebo") ~ "Placebo"
  ))

dplyr::select(metadata, refinebio_title, treatment_status)

str(metadata$treatment_status)

metadata <- metadata %>%
  dplyr::mutate(
    # Here we define the values our factor variable can have and their order.
    treatment_status = factor(treatment_status, levels = c("Stelara", "Placebo"))
  )

levels(metadata$treatment_status)

filteredExpressionDF <- expression_df %>%
  dplyr::filter(rowSums(.) >= 10)

geneMatrix <- round(filteredExpressionDF)
summary(geneMatrix)
any(geneMatrix < 0)

geneMatrix[geneMatrix < 0] <- 0


ddset <- DESeqDataSetFromMatrix(
  countData = geneMatrix,
  colData = metadata,
  design = ~treatment_status
)

deseq_object <- DESeq(ddset)

deseq_results <- results(deseq_object)
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)
head(deseq_results)

# this is of class DESeqResults -- we want a data frame
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))

head(deseq_df)

plotCounts(ddset, gene = "ENSG00000280670", intgroup = "treatment_status")

readr::write_tsv(
  deseq_df,
  file.path(
    results_dir,
    "GSE140684_diff_expr_results.tsv"
  )
)

# We'll assign this as `volcano_plot`
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)

# Print out plot here
volcano_plot

ggsave(
  plot = volcano_plot,
  file.path(plots_dir, "GSE140684_volcano_plot.png")
) # Replace with a plot name relevant to your data

# Run the DESeq2 analysis
dds <- DESeq(ddset)

# Get the results (for example, comparing two treatment groups)
res <- results(dds)

# Sort the results by adjusted p-value (or other ranking criteria like log2FoldChange)
res <- res[order(res$padj), ]

# Subset the top 50 genes
top50 <- head(res, 50)

# Convert to a data frame for easy export
top50_df <- as.data.frame(top50)

# View the top 50 genes table
head(top50_df)


# Save the full DESeq2 results to a CSV file in your results folder
write.csv(as.data.frame(res), "results/top50.csv", row.names = TRUE)
