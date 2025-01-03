##### Code to analyze diversity of microbiome dataset. Created by Katarzyna Morańska, 2024

#Install packages
# if (!requireNamespace("phyloseq", quietly = TRUE)) {
#   install.packages("BiocManager") 
#   BiocManager::install("phyloseq")  
# }
# 
#Install Biobase
# BiocManager::install("Biobase")

#CRAN Packages
# install.packages(c("tibble", "dplyr", "ggplot2"))

# Load libraries
library('phyloseq')
library('tibble')
library('ggplot2')
library('dplyr')
library('stats')
library('stats4')
library('Biobase')


###Files input
otu_mat <- read.delim('Otu_table_uterusvagina_t100_none_3_rare.tsv') #RAW abundance for alpha diversity
head(otu_mat)

samples <- read.delim("Meta_table_uterusvagina_t100.txt")
head(samples)

# Create phyloseq object
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



### Alpha diversity
otu_table(uteruseq) <- round(otu_table(uteruseq))

custom_colors <- c("uterus" = "#009E73", "vagina" = "#CC79A7", "H" = "#325a58")

plot_richness(uteruseq, measures=c("Chao1", "Shannon","Simpson"), x="Source", color="Diagnosis")+ geom_boxplot() + theme(axis.text.x = element_text(angle=0, vjust=0,6))+
  scale_color_manual(values = custom_colors)

richness = estimate_richness(uteruseq, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
richness

# Save the richness table as a TSV file
write.table(richness, file = "uterus_richness.tsv", sep = "\t", 
            row.names = TRUE, col.names = NA, quote = FALSE)


############## UTERUS AND VAGINA DIVERSITY_______________________________________
### Alpha diversity
otu_table(uteruseq) <- round(otu_table(uteruseq))

custom_colors <- c("uterus" = "#009E73", "vagina" = "#CC79A7", "H" = "#325a58")

plot_richness(uteruseq, measures=c("Chao1", "Shannon","Simpson"), x="Source", color="Source")+ geom_boxplot() + theme(axis.text.x = element_text(angle=0, vjust=0,6))+
  scale_color_manual(values = custom_colors)

plot_richness(uteruseq, measures=c("Chao1", "Shannon","Simpson"), x="Age", color="Source")+ geom_boxplot() + theme(axis.text.x = element_text(angle=0, vjust=0,6))+
  scale_color_manual(values = custom_colors)

plot_richness(uteruseq, measures=c("Chao1", "Shannon","Simpson"), x="BMI", color="Source")+ geom_boxplot() + theme(axis.text.x = element_text(angle=0, vjust=0,6))+
  scale_color_manual(values = custom_colors)

plot_richness(uteruseq, measures=c("Chao1", "Shannon","Simpson"), x="Menopause_status", color="Source")+ geom_boxplot() + theme(axis.text.x = element_text(angle=0, vjust=0,6))+
  scale_color_manual(values = custom_colors)


richness = estimate_richness(uteruseq, measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson"))
richness

# Save the richness table as a TSV file
write.table(richness, file = "uterusvagina_richness.tsv", sep = "\t", 
            row.names = TRUE, col.names = NA, quote = FALSE)









