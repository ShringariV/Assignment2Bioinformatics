# I want to visualise data using the bioinformatics data I have
import seaborn as sns

import matplotlib.pyplot as plt
import pandas as pd
# Assuming `data` is a dictionary with keys as sequences and values as their counts

file_path = 'oldData/metadata_GSE140684.tsv'
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

#Columns of the data
print(data.columns.tolist())

#for the columns of the metadata in ERP107715, generate in the age column, numbers between 5-12
