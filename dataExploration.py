# I want to visualise data using the bioinformatics data I have
import seaborn as sns

import matplotlib.pyplot as plt
import pandas as pd
# Assuming `data` is a dictionary with keys as sequences and values as their counts

metaDataFile = 'data/metadata_ERP107715.tsv'
metaData = pd.read_csv(metaDataFile, sep='\t')

#display first few rows
print(metaData.head())
print('\n')

#structure of dataframe
print("Data shape:")
print(metaData.shape)
print('\n')
print(metaData.info())
print('\n')

#Summary stats
print("Summary statistics:")
print(metaData.describe())
print('\n')

#Columns of the data
print(metaData.columns.tolist())

#for the columns of the metadata in ERP107715, generate in the age column, numbers between 5-12
dataFile = 'data/ERP107715.tsv'
data = pd.read_csv(dataFile, sep='\t')

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