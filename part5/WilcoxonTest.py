import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from scipy.stats import ranksums

df = pd.read_csv('../part3/ERP107715_diff_expr_results.tsv', sep='\t')

# print(df.columns)

upRegulated = df[df['log2FoldChange'] > 0]['log2FoldChange']
downRegulated = df[df['log2FoldChange'] < 0]['log2FoldChange']

stat, p_value = ranksums(upRegulated, downRegulated)

print(f'Wilcoxon rank-sum statistic: {stat}')
print(np.format_float_scientific(p_value, precision=100))

# Testing after getting zero p-value, I don't believe this is real.


print(f'Number of upregulated genes: {len(upRegulated)}')
print(f'Number of downregulated genes: {len(downRegulated)}')


plt.hist(upRegulated, bins=50, alpha=0.5, label='Upregulated')
plt.hist(downRegulated, bins=50, alpha=0.5, label='Downregulated')
plt.legend(loc='upper right')
plt.show()