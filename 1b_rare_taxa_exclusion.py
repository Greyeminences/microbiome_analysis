##### Code to proceed rare taxa exclusion on otu_table of microbiome dataset. Created by Katarzyna Morańska, 2024

import pandas as pd
import matplotlib.pyplot as plt

# Load the OTU table
### DATA IMPORT AND FILTERING
path_to_folder = r'paste path to the directory here'

# Path to files
source = 'uterus'               # 'uterus'/'uterusvagina'
threshold = '100'
normalization = 'normalizedCLR' # 'none'/'scaled'/'normalizedCLR'

otu_file = (f'{path_to_folder}\\Otu_table_{source}_t{threshold}_{normalization}.tsv')
otu_table = pd.read_csv(otu_file, sep='\t', index_col=0)

# Define thresholds for the minimum number of samples in which a taxon must appear
sample_presence_thresholds = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15,16,17,18,19,20]  # Adjust as needed
remaining_taxa_counts = []

# Filter taxa and count remaining taxa at each presence threshold
for sthreshold in sample_presence_thresholds:
    # Count the number of samples in which each taxon appears
    taxa_present_in_samples = (otu_table > 0).sum(axis=1)
    # Filter out taxa that are present in fewer than the threshold number of samples
    filtered_otu_table = otu_table.loc[taxa_present_in_samples >= sthreshold]
    remaining_taxa_counts.append(len(filtered_otu_table))

# Plot the number of taxa remaining at each presence threshold with labels
plt.figure(figsize=(10, 6))
plt.bar(sample_presence_thresholds, remaining_taxa_counts, color='purple')
plt.xlabel("Minimum number of samples with presence")
plt.ylabel("Number of taxa remaining")
plt.title(f'Number of taxa remaining after rare taxa exclusion in {source} t{threshold}')

# Add labels to each bar
for i, count in enumerate(remaining_taxa_counts):
    plt.text(sample_presence_thresholds[i], count, str(count), ha='center', va='bottom', fontsize=10, color='black')

plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()

# Save the filtered table at a specific presence threshold as TSV
selected_threshold = 3  # Edit threshold here
taxa_present_in_samples = (otu_table > 0).sum(axis=1)
filtered_otu_table = otu_table.loc[taxa_present_in_samples >= selected_threshold]
filtered_otu_table.to_csv(f'{path_to_folder}\\Otu_table_{source}_t{threshold}_{normalization}_{selected_threshold}_rare.tsv', sep='\t')
