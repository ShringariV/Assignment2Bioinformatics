# I want to visualise data using the bioinformatics data I have
import seaborn as sns

import matplotlib.pyplot as plt
import pandas as pd
# Assuming `data` is a dictionary with keys as sequences and values as their counts

file_path = 'BioinformaticsData/ERP107715/metadata_ERP107715.tsv'
data = pd.read_csv(file_path, sep='\t')

#display first few rows
print(data.head())
print('\n')

#structure of dataframe
print("Data shape:")
print(data.shape)
print('\n')
print(data.info())
print('\n')

#Summary stats
print("Summary statistics:")
print(data.describe())
print('\n')

#Missing values
print("Missing values:")
print(data.isnull().sum())
print('\n')

