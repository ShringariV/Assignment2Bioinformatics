results_dir <- file.path("results", "")

hgnFile <- file.path(results_dir, "GSE140684_Symbol.tsv")

expressionData <- read.csv(hgnFile, sep="\t", row.names = 1)

dim(expressionData)

non_numeric_values <- sapply(expressionData, class)
print(non_numeric_values)

expressionDataNumeric <- as.data.frame(lapply(expressionData, as.numeric))

sum(is.na(expressionDataNumeric))

expressionDataNumeric[is.na(expressionDataNumeric)] <- 0

logExpressionData <- log2(expressionDataNumeric+1)

head(logExpressionData)

geneMedians <- apply(logExpressionData, 1, median)

plot(density(geneMedians), 
     main = "Density of Per-Gene Median Expression (log2)", 
     xlab = "Log2 Median Expression", 
     ylab = "Density")
summary(geneMedians)
