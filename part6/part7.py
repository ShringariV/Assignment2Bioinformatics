import pandas as pd

# Load the joint enrichment results
input_file = "results/joint_enrichment_results.csv"
combined_df = pd.read_csv(input_file)

# Set a threshold for the number of methods to consider a term enriched in "most" methods
threshold = combined_df['methods_included'].max()  # Adjust if you want fewer methods

# Filter for terms enriched in most methods
top_terms_df = combined_df[combined_df['significant_count'] >= threshold - 1]  # Adjust as needed

# Sort by p-value to find the most significant terms
top_terms_df = top_terms_df.sort_values(by='pvalue').head(10)

# Save the top terms to a new CSV file
output_file = "results/top_enriched_terms.csv"
top_terms_df.to_csv(output_file, index=False)

print(f"Top enriched terms saved to {output_file}")