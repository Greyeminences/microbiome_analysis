##### Code to proceed otu_table normalization of microbiome dataset. Created by Rose Brouns and Katarzyna Morańska, 2024

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


### FUNCTIONS

def rename_cols(df):
    """
    Add 'sample_' to the columns names,
    except for the first column named '#Classification'
    """
    cols = df.columns.to_list()
    new_names = [cols[0]] + ['sample_' + col.split('.')[0] for col in cols[1:]]
    df.columns = new_names

    return(df)

def mean_overlap_cols(df):
    """
    Identify columns with the same name,
    and calculate the mean of the values.
    """
    overlapping_columns = df.columns[df.columns.duplicated()].unique()
    for col in overlapping_columns:
        cols_to_mean = df.loc[:, df.columns == col]
        df[col] = cols_to_mean.mean(axis=1)

    return (df)

def filter_otu_table_by_meta_table(fil_meta_table, otu_table):
    # Find samples Meta samples that are also in otu data
    samples = fil_meta_table['sample_name'].to_list()
    samples_fil = [col for col in samples if col in otu_table.columns]
    samples_fil = ['#Classification'] + samples_fil

    # Make new filtered OTU table
    fil_otu_table = otu_table[samples_fil]

    return(fil_otu_table)

def filter_meta_table_by_otu_table(meta_table, otu_table):
    samples = list(otu_table.columns.values)
    samples = list(filter(lambda name: name != "#Classification", samples))
    print(samples)

    return meta_table[meta_table['sample_name'].isin(samples)]

def filter_threshold_reads(df_num, threshold):
    # Calculate sum of values for each column
    col_sums = df_num.sum()

    # Filter columns where sum is greater than 800
    cols_over_threshold = col_sums[col_sums > threshold].index.tolist()

    return(cols_over_threshold)

def make_scaled_abundance(df_num):
    sum_reads = df_num.sum()
    df_relative_num = df_num.div(sum_reads, axis=1)

    return(df_relative_num)

def make_clr(df_pseudo):
    # Calculate the geometric mean for each sample (column-wise)
    geometric_means = np.exp(np.log(df_pseudo).mean(axis=0))

    # Calculate CLR for each feature in each sample
    df_clr = np.log(df_pseudo.div(geometric_means, axis=1))

    return(df_clr)


# input
path_to_folder = r'paste path to the directory here'
source = 'uterus'   # 'uterus'/'uterusvagina'

otu_file = (f'{path_to_folder}\\Otu_abundance_{source}.tsv')    #Edit name of the raw otu_file here
meta_file = (f'{path_to_folder}\\META_{source}.txt')            #Edit name of the meta_file here

# Read in tables
otu_table = pd.read_csv(otu_file, sep='\t')
meta_table = pd.read_csv(meta_file, sep='\t')

# # Filter meta_table by variable here if needed
# filtered_metatable = meta_table[meta_table['Diagnosis'] == 'H']
# meta_table = filtered_metatable.copy()

# # Add sample_name in meta_table to match otu_table names
# meta_table['sample_name'] = 'sample_' + meta_table['#NAME'].astype(str)

# Filter otu table by available samples in metadata 
renamed_otu = rename_cols(otu_table)

### STEP 1: Combine duplicates/triplicates

# Calculate mean of the overlapping columns
df_sum_cols = mean_overlap_cols(otu_table)

# Remove dupliacted cols
df_final = df_sum_cols.loc[:, ~df_sum_cols.columns.duplicated()]

# Filter otu_table by available samples in metadata 
fil_otu_table = filter_otu_table_by_meta_table(meta_table, df_final)
print(meta_table.shape[0])

# Filter meta_table by otu_table again
fil_meta_table = filter_meta_table_by_otu_table(meta_table, fil_otu_table)


### STEP 2: Normalization

# Get only numbers
df_num = fil_otu_table.iloc[:, 1:]
df_num_OG = df_num.copy()

# Filter on threshold
threshold =  100
samples_over_thr = filter_threshold_reads(df_num, threshold)
df_num = df_num[samples_over_thr]
df_num_final = pd.concat([df_final.iloc[:,0], df_num], axis=1)
df_num_final.fillna(0, inplace=True)

# Make scaled abundance
df_relative_num = make_scaled_abundance(df_num)
df_rela_final = pd.concat([df_final.iloc[:,0], df_relative_num], axis=1)
df_rela_final.fillna(0, inplace=True)

# Make CLR
# Make pseudocount for zeros
df_pseudo = df_num.replace(0, 1e-6)
df_clr = make_clr(df_pseudo)

# Add taxa
df_clr_final = pd.concat([df_final.iloc[:,0], df_clr], axis=1)

# Filter all meta table once again
fil_meta_table2 = filter_meta_table_by_otu_table(fil_meta_table, df_clr_final)


# Save to file
fil_otu_table.to_csv(f'{path_to_folder}\Otu_abundance_{source}.tsv', sep='\t', index=False)
fil_meta_table.to_csv(f'{path_to_folder}\Meta_abundance_{source}.txt', sep='\t', index=False)
df_num_final.to_csv(f'{path_to_folder}\Otu_table_{source}_t{threshold}_none.tsv', sep='\t', index=False)
df_rela_final.to_csv(f'{path_to_folder}\Otu_table_{source}_t{threshold}_scaled.tsv', sep='\t', index=False)
df_clr_final.to_csv(f'{path_to_folder}\Otu_table_{source}_t{threshold}_normalizedCLR.tsv', sep='\t', index=False)
fil_meta_table2.to_csv(f'{path_to_folder}\Meta_table_{source}_t{threshold}.txt', sep='\t', index=False)

# Variables' names match explanations
# fil_meta_table according to fil_otu_table
# fil_meta_table2 according to df_num_final (treshold), normalized data df_rela_final (scaled) and df_clr_final (CLR)

#  Files' names match explanations
# Meta_abundance[...] according to Otu_abundance[...]
# Meta_table[...] according to Otu_table[...]