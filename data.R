# Define the file paths
data_dir <- file.path("data")
data_file <- file.path(data_dir, "ERP107715.tsv")
metadata_file <- file.path(data_dir, "metadata_ERP107715.tsv")

# Check for the existence of files and import them
if (!file.exists(data_file)) {
  stop(paste("Data file does not exist:", data_file))
} else {
  data <- read.delim(data_file, header = TRUE, sep = "\t")
}

if (!file.exists(metadata_file)) {
  stop(paste("Metadata file does not exist:", metadata_file))
} else {
  metadata <- read.delim(metadata_file, header = TRUE, sep = "\t")
}

# Display the first few rows of the imported data
head(data)
head(metadata)