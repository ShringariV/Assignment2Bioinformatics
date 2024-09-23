import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import ranksums

#Gene-Ontology
df = pd.read_csv('../part3/ERP107715_diff_expr_results.tsv', sep='\t')

upRegulated = df[df['log2FoldChange'] > 0]['log2FoldChange']
downRegulated = df[df['log2FoldChange'] < 0]['log2FoldChange']

stat, p_value = ranksums(upRegulated, downRegulated)

# save stat and p_value to a csv
pd.DataFrame([{"statistic": stat, "p_value": np.format_float_scientific(p_value, precision=100)}]).to_csv(
    "enrichment results/wilcoxon_results.csv", index=False)

# Testing after getting zero p-value, I don't believe this is real.

print(f'Number of upregulated genes: {len(upRegulated)}')
print(f'Number of downregulated genes: {len(downRegulated)}')

plt.hist(upRegulated, bins=50, alpha=0.5, label='Upregulated')
plt.hist(downRegulated, bins=50, alpha=0.5, label='Downregulated')
plt.legend(loc='upper right')

# save plot to a file
plt.savefig("genes_plot.png")

plt.show()