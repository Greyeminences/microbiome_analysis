import pandas as pd
from sklearn.manifold import MDS
import matplotlib.pyplot as plt

# Load data 
meta_df = pd.read_csv('Meta_table_uterus_t100.txt', sep='\t')
dissimilarity_matrix = pd.read_csv("uterus_unifrac.tsv", sep=",", index_col=0)

# Perform MDS
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
embedding_weighted = mds.fit_transform(dissimilarity_matrix)

# Convert embeddings to DataFrames
embedding_df_weighted = pd.DataFrame(embedding_weighted, index=dissimilarity_matrix.columns, columns=["MDS1", "MDS2"])

# Merge with metadata
embedding_df_weighted = embedding_df_weighted.merge(meta_df[['sample_name', 'Microbiome_profile', 'NAME']], 
                                                    left_index=True, right_on="sample_name")

# Define colors for each Microbiome_profile category
color_map = {1: "darkblue", 2: "lightblue", 3: "green", 4: "pink"}
embedding_df_weighted["color"] = embedding_df_weighted["Microbiome_profile"].map(color_map)

# Plotting for weighted
fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(embedding_df_weighted["MDS1"], embedding_df_weighted["MDS2"], 
           c=embedding_df_weighted["color"], s=150, edgecolors='black')

for row in embedding_df_weighted.itertuples(index=False):
    x, y, label_name = row.MDS1, row.MDS2, row.NAME
    ax.text(x + 0.01, y + 0.02, label_name, fontsize=10, ha="center", va="bottom")

ax.set_xlabel("MDS 1")
ax.set_ylabel("MDS 2")
ax.set_title("Weighted MDS Projection UTERUS")
ax.set_aspect('equal')

# Legend manually
handles = []
for profile, color in color_map.items():
    # Create a dummy line for each profile to use in the legend
    handle = plt.Line2D([0], [0], marker='o', color='w', label=f"Type {profile}", 
                         markerfacecolor=color, markersize=10)
    handles.append(handle)
ax.legend(handles=handles, title="Microbiome type", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.show()



#####uterus vagina##########################
# Load data 
meta_df = pd.read_csv('Meta_table_uterusvagina_t100.txt', sep='\t')
dissimilarity_matrix = pd.read_csv("uterusvagina_unifrac.tsv", sep=",", index_col=0)

# Perform MDS
mds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
embedding_weighted = mds.fit_transform(dissimilarity_matrix)

# Convert embeddings to DataFrames
embedding_df_weighted = pd.DataFrame(embedding_weighted, index=dissimilarity_matrix.columns, columns=["MDS1", "MDS2"])

# Merge with metadata
embedding_df_weighted = embedding_df_weighted.merge(meta_df[['sample_name', 'Microbiome_profile', 'NAME', 'Patient_code', 'Source']], 
                                                    left_index=True, right_on="sample_name")

# Define colors for each Microbiome_profile category
color_map = {1: "darkblue", 2: "lightblue", 3: "green", 4: "pink"}
embedding_df_weighted["color"] = embedding_df_weighted["Microbiome_profile"].map(color_map)

# Plotting for weighted
fig, ax = plt.subplots(figsize=(10, 10))

# Uteruses
uteruses = embedding_df_weighted[embedding_df_weighted["Source"] == 'uterus']
ax.scatter(uteruses["MDS1"], uteruses["MDS2"], 
           c=uteruses["color"], s=150, marker='o', edgecolors='black')

# Vaginas
vaginas = embedding_df_weighted[embedding_df_weighted["Source"] == 'vagina']
ax.scatter(
    vaginas["MDS1"], 
    vaginas["MDS2"], 
    c=vaginas["color"], s=150, marker='v', edgecolors='black')

# Annotate points with "NAME"
for row in embedding_df_weighted.itertuples(index=False):
    x, y, label_name = row.MDS1, row.MDS2, row.NAME
    ax.text(x + 0.01, y + 0.02, label_name, fontsize=10, ha="center", va="bottom")

# Group the data by 'Patient_code' to find pairs of points
grouped_weighted = embedding_df_weighted.groupby('Patient_code')

# Draw lines between points with the same Patient_code and add the Patient_code label on the line
for patient, group in grouped_weighted:
    if len(group) > 1:
        for i in range(len(group) - 1):
            p1 = group.iloc[i]
            p2 = group.iloc[i + 1]
            
            # Draw line between the two points
            ax.plot(
                [p1['MDS1'], p2['MDS1']], 
                [p1['MDS2'], p2['MDS2']], 
                color='gray', 
                linestyle='--', 
                linewidth=1
            )

ax.set_xlabel("MDS 1")
ax.set_ylabel("MDS 2")
ax.set_title("Weighted MDS Projection UTERUS & VAGINA")
ax.set_aspect('equal')

# Add legend with labels and color mapping
ax.legend(title="Patient Code", bbox_to_anchor=(1.05, 1), loc='upper left')

# Legend 1: Microbiome profile (color only)
profile_handles = []
for profile, color in color_map.items():
    handle = plt.Line2D(
        [0], [0], marker='o', color='w', label=f"Type {profile}",
        markerfacecolor=color, markeredgecolor='black', markersize=10
    )
    profile_handles.append(handle)

legend1 = ax.legend(handles=profile_handles, title="Microbiome type", loc='upper left', bbox_to_anchor=(1.05, 1))
ax.add_artist(legend1)  # Add first legend manually so it's not overwritten

# Legend 2: Source (shape only)
source_handles = [
    plt.Line2D([0], [0], marker='o', color='w', label='Uterus',
               markerfacecolor='gray', markeredgecolor='black', markersize=10),
    plt.Line2D([0], [0], marker='v', color='w', label='Vagina',
               markerfacecolor='gray', markeredgecolor='black', markersize=10),
]

ax.legend(handles=source_handles, title="Sample source", loc='upper left', bbox_to_anchor=(1.05, 0.8))


plt.show()