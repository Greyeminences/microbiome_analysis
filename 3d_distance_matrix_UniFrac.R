###UNIFRAC

# Load libraries
library(ape)
library(phyloseq)
library(dplyr)

### Load the tree
tree <- read.tree("tree_file.tree")
print(tree)
summary(tree)
#plot(tree)
tree_taxa <- tree$tip.label
head(tree_taxa)

### Load otu file and create phyloseq object
otu <- read.delim('uterus_otu_file_withtreecorrespondinglabels.tsv')
head(otu)

merged_otu <- otu %>%
  group_by(Classification) %>%  
  summarise(across(everything(), sum, na.rm = TRUE))

merged_otu <- merged_otu %>%
  tibble::column_to_rownames("Classification") 

merged_otu <- as.matrix(merged_otu)

physeq <- phyloseq(otu_table(merged_otu, taxa_are_rows = TRUE), phy_tree(tree))
physeq
physeq@phy_tree

# Calculate weighted UniFrac distance matrix
weighted_uniFrac_dist <- UniFrac(physeq, weighted = TRUE)
weighted_uniFrac_dist

# Convert UniFrac distance matrix to a data frame
weighted_uniFrac_df <- as.data.frame(as.matrix(weighted_uniFrac_dist))

# View the first few rows of the data frame
head(weighted_uniFrac_df)

# Save the data frame as a CSV file
write.csv(weighted_uniFrac_df, "uterus_unifrac.tsv", row.names = TRUE)