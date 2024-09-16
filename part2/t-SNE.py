import pandas as pd
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns

expressionData = pd.read_csv('../data/GSE140684_Symbol.tsv', sep='\t', index_col=0)

metaData = pd.read_csv('../data/metadata_GSE140684.tsv', sep='\t', index_col=0)

commonSamples = expressionData.columns.intersection(metaData.index)
expressionData = expressionData[commonSamples]
metaData = metaData.loc[commonSamples]

expressionMatrix = expressionData.values

expressionMatrix = np.log1p(expressionMatrix)

tsne = TSNE(n_components=2, random_state=42)
tsne_results = tsne.fit_transform(expressionMatrix.T)

metaData['tsne_2d-one'] = tsne_results[:,0]
metaData['tsne_2d-two'] = tsne_results[:,1]

plt.figure(figsize=(10,8))
sns.scatterplot(
    x="tsne_2d-one", y="tsne_2d-two",
    hue="refinebio_treatment",
    palette=sns.color_palette("viridis", as_cmap=False),
    data=metaData,
    legend="full",
    alpha=0.8
)

plt.title('t-SNE plot of Gene Expression Data')

plt.savefig('tsne_plot.png', dpi=300)

plt.show()