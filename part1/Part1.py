import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load the data from the TSV file into a pandas DataFrame
file_path = "../data/ERP107715_GeneNames.tsv"
df = pd.read_csv(file_path, sep='\t')

#print(df.head())

expressionData = df.iloc[:, 3:]

#log transformation
logExpressionData = np.log2(expressionData + 1)

#median expression per gene
medianExpression = logExpressionData.median(axis=1)

#density plot
sns.set(style="whitegrid")
plt.figure(figsize=(8, 6))
sns.kdeplot(medianExpression,fill=True,color='blue',alpha=0.5)
plt.title("Density of Per-Gene Median Expression (Log-Scaled)")
plt.xlabel("Median Expression (log2 scale)")
plt.ylabel("Density")
plt.show()