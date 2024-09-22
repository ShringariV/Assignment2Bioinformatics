import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

expressionData = pd.read_csv("../data/ERP107715.tsv", sep='\t')
print(expressionData.head())

if expressionData is not None and 'Gene' in expressionData.columns:
    expressionData = expressionData.set_index('Gene')  # Reassign to get the modified DataFrame
else:
    raise ValueError("The expression data could not be loaded or does not contain a 'Gene' column.")

metadata = pd.read_csv("../data/metadata_ERP107715.tsv", sep='\t')

conditions = metadata['refinebio_sex']

pValues = []

logFoldChanges = []

for Gene in expressionData.index:
    print("Processing gene " + Gene + "")
    print("\n")

    data = pd.DataFrame({
        'expression': expressionData.loc[Gene].values,
        'condition': conditions
    })
    conditionOneExpression = data.loc[data['condition'] == data['condition'].unique()[0], 'expression']
    conditionTwoExpression = data.loc[data['condition'] == data['condition'].unique()[1], 'expression']

    p = mannwhitneyu(conditionOneExpression, conditionTwoExpression)[1]
    logFoldChange = np.log2((conditionTwoExpression.mean() +1) / (conditionOneExpression.mean() + 1))
    pValues.append(p)
    logFoldChanges.append(logFoldChange)

fdr = multipletests(pValues, method='fdr_bh')[1]

results = pd.DataFrame({
    'gene': expressionData.index,
    'pValue': pValues,
    'fdr': fdr,
    'logFoldChange': logFoldChanges
})

print(results)

#save this into a file
results.to_csv("../data/mannWhitneyUTest_results.csv", index=False)