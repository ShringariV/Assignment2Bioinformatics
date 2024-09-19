import pandas as pd

# Load the metadata file
metadata_df = pd.read_csv('data/metadata_ERP107715.tsv', sep='\t')

# Update 'refinebio_title' with the combined text from 'refinebio_title' and 'refinebio_sex'
metadata_df['refinebio_title'] = metadata_df['refinebio_title'] + ' - ' + metadata_df['refinebio_sex']

# Save the updated DataFrame to a new TSV file
metadata_df.to_csv('updated_metadata_ERP107715.tsv', sep='\t', index=False)
