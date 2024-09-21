# Volcano Plot Summary

- **Title**: EnhancedVolcano
- **X-Axis**: Log2 fold change
- **Y-Axis**: -Log10 P
- **Data Points**:
  - Blue: Not Significant (NS)
  - Green: Significant by Log2 Fold Change
  - Red: Significant by p-value
  - Purple: Significant by both p-value and Log2 Fold Change
- **Significance Thresholds**:
  - Vertical Lines at Log2 fold change values of -2 and 2.
  - Horizontal Line at a p-value threshold of approximately \(10^{-50}\).
- **Labeled Points**:
  - Top Right: ENSG00000129824 (Significant by both criteria)
  - Bottom Left: ENSG00000229807 (Not significant)
- **Total Variables**: 29708


# Gene Expression Data Summary

- **Columns**:
  - **Gene**: Identifier for each gene.
  - **log2FoldChange**: Log2-transformed fold change in gene expression.
  - **lfcSE**: Standard error of the log2 fold change.
  - **stat**: Test statistic.
  - **pvalue**: P-value indicating the significance of the gene expression change.
  - **padj**: Adjusted p-value for multiple testing correction.
