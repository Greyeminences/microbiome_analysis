##### Code to proceed dimension reduction (PCA, tSNE) of microbiome dataset. Created by Rose Brouns and Katarzyna Morańska, 2024

import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import datetime
import matplotlib.pyplot as plt
import math

from pathlib import Path


### FUNCTIONS

def filter_meta_table_excluded_samples(meta_table):
    # Drop excluded samples from analysis
    fil_meta_table = meta_table[~meta_table.apply(lambda row: row.astype(str).str.contains('excluded').any(), axis=1)]

    return(fil_meta_table)

def filter_otu_table_by_meta_table(fil_meta_table, otu_table):
    # Find samples Meta samples that are also in otu data
    samples = fil_meta_table['sample_name'].to_list()
    samples_fil = [col for col in samples if col in otu_table.columns]
    samples_fil = ['#Classification'] + samples_fil

    # Make new filtered OTU table
    fil_otu_table = otu_table[samples_fil]

    return(fil_otu_table)

def filter_otu_table_by_meta_table_v2(fil_meta_table, otu_table):
    # Find samples Meta samples that are also in otu data
    samples = fil_meta_table['#NAME'].to_list()
    samples_fil = ["# "+str(col) for col in samples if "# "+str(col) in otu_table.columns]
    samples_fil = ['#Classification'] + samples_fil

    # Make new filtered OTU table
    fil_otu_table = otu_table[samples_fil]

    return(fil_otu_table)

def add_water_control_col(df):
    df['water_control'] = np.where(df['sample_name'].str.contains('water', case=False, na=False), 'yes','no')
    return(df)

def plot_pca(pca_df, hue_opt, color_palette=None):
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=pca_df['PC1'], y=pca_df['PC2'], alpha=0.5, hue=pca_df[hue_opt], palette=color_palette)

    plt.legend(title='Color Column', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title(f'PCA Plot {hue_opt}')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')    
    plt.grid(True)
    plt.show()

def plot_tsne(df_tsne, hue_opt, title_opt='', color_palette=None):
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=df_tsne['t-sne_1'], y=df_tsne['t-sne_2'], hue=df_tsne[hue_opt], palette=color_palette, alpha=0.5)

    plt.legend(title='Color Column', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title(f't-SNE Plot {hue_opt} {title_opt}')
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.grid(True)
    plt.show()

def rename_cols(df):
    """
    Add 'sample_' to the columns names,
    except for the first column named '#Classification'
    """
    cols = df.columns.to_list()
    new_names = [cols[0]] + ['sample_' + col.split('.')[0] for col in cols[1:]]
    df.columns = new_names

    return(df)


### DATA IMPORT AND FILTERING
path_to_folder = r'#paste path to the directory here'

# Path to files
source = 'uterus'                   # 'uterus'/'uterusvagina'
threshold = '100'
normalization = 'normalizedCLR'     # 'none'/'scaled'/'normalizedCLR'
rare_excluded = '3'

otu_file = (f'{path_to_folder}\\Otu_table_{source}_t{threshold}_{normalization}_{rare_excluded}_rare.tsv')
meta_file = (f'{path_to_folder}\\Meta_table_{source}_t{threshold}.txt')

# Read in tables
otu_table = pd.read_csv(otu_file, sep='\t')
meta_table = pd.read_csv(meta_file, sep='\t')

# Set colors
custom_palette = {
    'E': '#ADADAD',
    'O': '#ADADAD',
    'H': "#325a58",
    'OC': "#D55E00", 
    'EC': "#0072B2",  
    'after': "#D81B60",  
    'before': "#004D40",
    'Normal' : "#027D5C",
    'Overweight' : "#E2BC00",
    'Obesity' : "#F12929",  
    "uterus" : "#009E73", 
    "vagina" : "#CC79A7" 
}

colors = custom_palette


### PCA

# Set PCA parameters
pca = PCA(n_components=2, random_state=42)

# All data PCA
all_matrix = otu_table.iloc[:,1:].T
pca_result = pca.fit_transform(all_matrix)
pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
pca_df['sample_name'] = otu_table.columns[1:]

# Add metadata to PCA df
pca_df_meta = pca_df.merge(meta_table, how='left', on='sample_name')

# Explained variance ratio
explained_variance = pca.explained_variance_ratio_
explained_variance_ratio = pca.explained_variance_ratio_
print("Explained variance ratio:", explained_variance_ratio)
cumulative_explained_variance = explained_variance_ratio.cumsum()
print("Cumulative explained variance:", cumulative_explained_variance)


### PCA in a loop by variable

for variable in meta_table.columns [[4,6,7,8,9,15]]:   
    # Filter the dataframe to exclude missing values for the current variable
    filtered_pca_df_meta = pca_df_meta[pca_df_meta[variable].notna()]
      
    fig = plt.figure(figsize=(16, 12))
    fig.subplots_adjust(hspace=0.2)
    color_palette = colors if variable in ["Source", "Menopause_status", "BMI", "Diagnosis"] else None
    
    sns.scatterplot(x=filtered_pca_df_meta['PC1'], y=filtered_pca_df_meta['PC2'], hue=filtered_pca_df_meta[variable], palette=color_palette, s=300, alpha=1)
    
    # Add labels to each point
    for i in list(filtered_pca_df_meta.index):

        # Fisrt label: Microbiome_profile
        plt.text(
            filtered_pca_df_meta['PC1'][i], 
            filtered_pca_df_meta['PC2'][i], 
            filtered_pca_df_meta['Microbiome_profile'][i], 
            fontsize=10, 
            ha='center', 
            va='center',
            color='white',
            fontweight='bold'
        )

        # Second label: source_name
        plt.text(
            filtered_pca_df_meta['PC1'][i] + 2.5, 
            filtered_pca_df_meta['PC2'][i],
            filtered_pca_df_meta['source_name'][i], 
            fontsize=10, 
            ha='left', 
            va='center', 
            color='black', 
            fontweight='normal'
        )
    
    plt.legend(title='Color Column', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title(f'PCA Plot {variable}')
    plt.xlabel(f'PC 1 ({explained_variance[0]:.2%})')
    plt.ylabel(f'PC 2 ({explained_variance[1]:.2%})')
    plt.grid(True)
    plt.show()

    # Generate a timestamp
    timestamp = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
    plt.tight_layout()
    
    # Define the path to the output folder
    output_folder = Path(f'{path_to_folder}/Output/plots')

    # Create the directory if it doesn't exist
    output_folder.mkdir(parents=True, exist_ok=True)

    # Save the plot
    fig.savefig(output_folder / f'PCA_{source}_{variable}_t{threshold}_{normalization}_{rare_excluded}_rare_{timestamp}.png', bbox_inches='tight')
    fig.savefig(output_folder / f'PCA_{source}_{variable}_t{threshold}_{normalization}_{rare_excluded}_rare_{timestamp}_svg.svg', bbox_inches='tight')




####### PCA FOR UTERUS AND VAGINA WITH DASHED LINE BETWEEN SOURCES____________________________________________

### PCA in a loop by variable
for variable in meta_table.columns [[4,6,7,8,9,15]]:   
    # Filter the dataframe to exclude missing values for the current variable
    filtered_pca_df_meta = pca_df_meta[pca_df_meta[variable].notna()]
      
    fig = plt.figure(figsize=(16, 12))
    fig.subplots_adjust(hspace=0.2)
    color_palette = colors if variable in ["Source", "Menopause_status", "BMI", "Diagnosis"] else None
    
    sns.scatterplot(x=filtered_pca_df_meta['PC1'], y=filtered_pca_df_meta['PC2'], hue=filtered_pca_df_meta[variable], palette=color_palette, s=300, alpha=1)
    
    # Add labels to each point
    for i in list(filtered_pca_df_meta.index):

        # Fisrt label: Microbiome_profile
        plt.text(
            filtered_pca_df_meta['PC1'][i], 
            filtered_pca_df_meta['PC2'][i], 
            filtered_pca_df_meta['Microbiome_profile'][i], 
            fontsize=10, 
            ha='center', 
            va='center',
            color='white',
            fontweight='bold'
        )

    # Group the data by 'Patient_code' to find pairs of points
    grouped = filtered_pca_df_meta.groupby('Patient_code')
    
    # Draw lines between points with the same Patient_code and add the Patient_code label on the line
    for patient, group in grouped:
        if len(group) > 1:
            for i in range(len(group) - 1):
                p1 = group.iloc[i]
                p2 = group.iloc[i + 1]
                
                # Draw line between the two points
                plt.plot(
                    [p1['PC1'], p2['PC1']], 
                    [p1['PC2'], p2['PC2']], 
                    color='gray', 
                    linestyle='--', 
                    linewidth=1
                )
                
                # Calculate midpoint of the line
                midpoint_x = (p1['PC1'] + p2['PC1']) / 2
                midpoint_y = (p1['PC2'] + p2['PC2']) / 2
                
                # Calculate the angle of the line (in radians) using np.arctan2
                delta_y = p2['PC2'] - p1['PC2']
                delta_x = p2['PC1'] - p1['PC1']
                angle = np.arctan2(delta_y, delta_x)  # This gives the angle of the line in radians
                angle_deg = np.degrees(angle)  # Convert to degrees

                angle_threshold = 90

                # Determine the rotation condition based on the threshold
                if abs(angle_deg) > angle_threshold:
                    rotation_angle = angle_deg  # Rotate the label if the angle exceeds the threshold
                else:
                    rotation_angle = 0  # Keep the label horizontal if below the threshold

                # Add Patient_code label at the midpoint of the line, rotated along the line
                plt.text(
                    midpoint_x, 
                    midpoint_y, 
                    patient,   # Patient_code label
                    fontsize=10, 
                    ha='center', 
                    va='center', 
                    color='black',  # Choose a contrasting color for visibility
                    fontweight='normal',
                    rotation=np.degrees(angle)  # Convert radians to degrees for rotation
                )

    plt.legend(title='Color Column', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title(f'PCA Plot {variable}')
    plt.xlabel(f'PC 1 ({explained_variance[0]:.2%})')
    plt.ylabel(f'PC 2 ({explained_variance[1]:.2%})')
    plt.grid(True)
    plt.show()

    # Generate a timestamp
    timestamp = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
    plt.tight_layout()
    
    # Define the path to the output folder
    output_folder = Path(f'{path_to_folder}/Output/plots')

    # Create the directory if it doesn't exist
    output_folder.mkdir(parents=True, exist_ok=True)

    # Save the plot
    fig.savefig(output_folder / f'PCA_{source}_{variable}_t{threshold}_{normalization}_{rare_excluded}_rare_{timestamp}.png', bbox_inches='tight')
    fig.savefig(output_folder / f'PCA_{source}_{variable}_t{threshold}_{normalization}_{rare_excluded}_rare_{timestamp}_svg.svg', bbox_inches='tight')

#####___________________________________________________________________________________________________________________


### t-SNE
tsne = TSNE(n_components=2, perplexity = 21, random_state=42)

# Whole dataset
X_tsne = tsne.fit_transform(all_matrix)
df_tsne = pd.DataFrame(X_tsne, columns=['t-sne_1', 't-sne_2'])
df_tsne['sample_name'] = otu_table.columns[1:]

# Add metadata
df_tsne_meta = df_tsne.merge(meta_table, how='left', on='sample_name')


## tSNE in a loop per variable
for perplexity_value in [5, math.floor(len(all_matrix)/2), len(all_matrix)-1]:
    
    for variable in meta_table.columns [[4,6,7,8,9,15]]: 

        tsne = TSNE(n_components=2, perplexity=perplexity_value, random_state=42)
        X_tsne = tsne.fit_transform(all_matrix)
    
        df_tsne = pd.DataFrame(X_tsne, columns=['t-sne_1', 't-sne_2'])
        df_tsne['sample_name'] = otu_table.columns[1:]
      
        # Add metadata
        df_tsne_meta = df_tsne.merge(meta_table, how='left', on='sample_name')

        # Filter samples without metadata for the current variable
        filtered_df_tsne_meta = df_tsne_meta[df_tsne_meta[variable].notna()]

        fig = plt.figure(figsize=(16, 14))
        fig.subplots_adjust(hspace=0.2)
        
        color_palette = colors if variable in ["Source", "Menopause_status", "BMI", "Diagnosis"] else None
    
        sns.scatterplot(x=filtered_df_tsne_meta['t-sne_1'], y=filtered_df_tsne_meta['t-sne_2'], hue=filtered_df_tsne_meta[variable], palette=color_palette, s=300, alpha=1)
        
        # Add labels to each point
        for i in list(filtered_df_tsne_meta.index):

            # Fisrt label: Microbiome_profile
            plt.text(
                filtered_df_tsne_meta['t-sne_1'][i], 
                filtered_df_tsne_meta['t-sne_2'][i], 
                filtered_df_tsne_meta['Microbiome_profile'][i], 
                fontsize=10, 
                ha='center', 
                va='center',
                color='white',
                fontweight='bold'
            )

            # Second label: source_name
            plt.text(
                filtered_df_tsne_meta['t-sne_1'][i] + 2.5, 
                filtered_df_tsne_meta['t-sne_2'][i],
                filtered_df_tsne_meta['source_name'][i], 
                fontsize=10, 
                ha='left', 
                va='center', 
                color='black', 
                fontweight='normal'
            )

        plt.legend(title='Color Column', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.title(f't-SNE Plot {variable}, perplexity {perplexity_value}')
        plt.xlabel('t-SNE 1')
        plt.ylabel('t_SNE 2')
        plt.grid(True)
        plt.show()

        # Generate a timestamp
        timestamp = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
        plt.tight_layout()

        # Define the path to the folder
        output_folder = Path(f'{path_to_folder}/Output/plots')

        # Create the directory if it doesn't exist
        output_folder.mkdir(parents=True, exist_ok=True)

        # Save the plot
        fig.savefig(output_folder / f'tSNE_{source}_{variable}_t{threshold}_{normalization}_{rare_excluded}_rare_{timestamp}.png', bbox_inches='tight')
        fig.savefig(output_folder / f'tSNE_{source}_{variable}_t{threshold}_{normalization}_{rare_excluded}_rare_{timestamp}_svg.svg', bbox_inches='tight')
