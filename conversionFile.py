import pandas as pd
from pybiomart import Dataset

df = pd.read_csv('part3/ERP107715_diff_expr_results.tsv', sep='\t')

nulls = df['threshold'].isnull()

df.loc[nulls, 'threshold'] = 'False'

df.to_csv('differential_expression_results.tsv', sep='\t', index=False)
