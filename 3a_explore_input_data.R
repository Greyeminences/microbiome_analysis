##### Code to explore raw input seq data. Created by Katarzyna Morańska, 2024

##Install packages
# # Install BiocManager if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# # Install phyloseq and Bioconductor dependencies
# BiocManager::install(c("phyloseq"))
# 
# # Install tidyverse (includes tibble, ggplot2, and dplyr)
# install.packages("tidyverse")

# Load libraries
library('tibble')
library('phyloseq')
library('ggplot2')
library('dplyr')


###Files input
otu_mat <- read.delim('Otu_table_uterus_t100_none_3_rare.tsv')
head(otu_mat)

samples <- read.delim("Meta_table_uterus_t100.txt")
head(samples)

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("X.Classification") 

samples <- samples %>% 
  tibble::column_to_rownames("sample_name") 

otu_mat <- as.matrix(otu_mat)

# Extract the rownames (SILVA long names) from the OTU table
silva_long_names <- rownames(otu_mat)

# Split the SILVA long names into separate taxonomic ranks by the '|' delimiter
taxonomy <- strsplit(silva_long_names, "\\|") %>%
  lapply(function(x) {
    # Remove the "k__", "p__", etc. and keep only the taxon names
    sapply(x, function(level) sub("^[a-z]__", "", level))  # This removes k__, p__, c__, etc.
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

# Assign column names to represent taxonomic ranks
colnames(taxonomy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:ncol(taxonomy)]

# Set rownames (corresponding to OTUs) to match the OTU table
rownames(taxonomy) <- rownames(otu_mat)

# Convert to matrix for phyloseq
tax_matrix <- as.matrix(taxonomy)

# Check the first few rows of the generated taxonomy table
head(tax_matrix)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
samples = sample_data(samples)

uteruseq <- phyloseq(OTU, TAX, samples)
uteruseq

# Make relative abundance if needed
#uteruseq_relative <- transform_sample_counts(uteruseq, function(x) x / sum(x))  


### plot nb of raw reads
plot_bar(uteruseq, y = "Abundance", fill = "Genus") +
  xlab("sample_id") +
  ylab("Raw abundance") +
  #geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  #scale_fill_manual(values = unique_colors_df_vector)+
  geom_hline(yintercept = 100, color = "red", linetype = "dashed", size = 1) + 
  theme(legend.position = "none")
