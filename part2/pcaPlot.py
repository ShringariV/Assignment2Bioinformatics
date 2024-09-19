import os

import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# Load your expression data and metadata
try:
    if os.path.exists('../data/ERP107715_Symbol.tsv'):
        expression_data = pd.read_csv('../data/ERP107715_Symbol.tsv', sep='\t')
    else:
        print("File ../data/ERP107715_Symbol.tsv does not exist")

    if os.path.exists('../oldData/metadataERP107715Filtered.tsv'):
        metadata = pd.read_csv('../oldData/metadataERP107715Filtered.tsv', sep='\t')
    else:
        print("File ../data/metadataERP107715Filtered.tsv does not exist")

except pd.errors.ParserError:
    print("Error: File is not in the correct format")
except PermissionError:
    print("Error: No permissions to read the file")

# Ensure consistent filtering between expression data and metadata
common_columns = [col for col in expression_data.columns[3:] if col in metadata['Comment[ENA_RUN]'].values]
expression_data_filtered = expression_data[common_columns]
metadata_filtered = metadata[metadata['Comment[ENA_RUN]'].isin(common_columns)]

print(common_columns)
print(expression_data_filtered)
print(metadata_filtered)

X = expression_data_filtered.iloc[:, 3:].values
print(X)

xStandardised = StandardScaler().fit_transform(X.T)

print(xStandardised)

pca = PCA(n_components=2)
principalComponents = pca.fit_transform(xStandardised)

print(principalComponents)
print(pca)

# Create a DataFrame for the principal components
pca_df = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])

# Check the length of common_columns and pca_df.index
print(f"Length of common_columns: {len(common_columns)}")
print(f"Length of PCA output: {len(pca_df)}")

# Ensure the lengths match
if len(common_columns) == len(metadata_filtered):
    print("Columns and metadata lengths match!")
else:
    print(f"Mismatch: {len(common_columns)} expression columns vs {len(metadata_filtered)} metadata entries.")
pca_df = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])