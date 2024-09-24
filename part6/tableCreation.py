import os
import pandas as pd

# Directory where the enrichment result CSV files are located
enrichment_results_dir = "../part5/enrichment results"

# Create an empty list to hold the dataframes from each CSV file
dfs = []

# Loop through each CSV file in the directory
for file_name in os.listdir(enrichment_results_dir):
    if file_name.endswith(".csv"):
        file_path = os.path.join(enrichment_results_dir, file_name)
        df = pd.read_csv(file_path)

        # Add a column to indicate which method this data is from (based on file name)
        df['method'] = file_name.split('.')[0]

        # Append to the list of dataframes
        dfs.append(df)

# Combine all dataframes into one large dataframe
combined_df = pd.concat(dfs)

# Group by the gene set/term and calculate how many methods found each term significantly enriched
# Assuming that 'padj' is the adjusted p-value column for significance threshold (e.g., < 0.05)
combined_df['significant'] = combined_df['padj'] < 0.05

# Group by the gene set/term to get unique terms and their significance counts
summary_df = combined_df.groupby('id').agg(
    methods_included=('method', 'nunique'),  # Number of methods that included the term
    significant_count=('significant', 'sum')  # Number of methods that found the term significant
).reset_index()

# Merge with all the relevant statistics from each method
final_df = pd.merge(combined_df, summary_df, on='id', how='left')

# Save the final table to the results folder
output_file = os.path.join("results", "joint_enrichment_results.csv")
final_df.to_csv(output_file, index=False)

print(f"Joint table of enriched terms saved to {output_file}")
