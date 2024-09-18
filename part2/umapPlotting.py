import pandas as pd
import numpy as np
import umap
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
import umapPlotting

expressionData = pd.read_csv('../oldData/GSE140684_Symbol.tsv', sep='\t', index_col=0)

metaData = pd.read_csv('../oldData/metadata_GSE140684.tsv', sep='\t', index_col=0)

commonSamples = expressionData.columns.intersection(metaData.index)
expressionData = expressionData[commonSamples]
metaData = metaData.loc[commonSamples]

expressionMatrix = expressionData.values

expressionMatrix = np.log1p(expressionMatrix)

umap_model = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
umap_results = umap_model.fit_transform(expressionMatrix.T)

metaData['umap-2d-one'] = umap_results[:, 0]
metaData['umap-2d-two'] = umap_results[:, 1]

# Create and save UMAP plot
plt.figure(figsize=(10, 8))
sns.scatterplot(
    x="umap-2d-one", y="umap-2d-two",
    hue="refinebio_treatment",
    palette=sns.color_palette("viridis", as_cmap=False),
    data=metaData,
    legend="full",
    alpha=0.8
)

plt.title('UMAP plot of Gene Expression Data')

# Save the plot to a file
plt.savefig('umap_plot.png', dpi=300)  # Save as PNG with 300 DPI resolution

# Optionally, display the plot
plt.show()