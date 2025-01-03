##### Code to extract from microbiome dataset top n taxa by abundance. Created by Katarzyna Morańska, 2024

#Install packages
# if (!requireNamespace("phyloseq", quietly = TRUE)) {
#   install.packages("BiocManager")  # Install BiocManager if needed
#   BiocManager::install("phyloseq")  # Install phyloseq
# }
# 
# if (!requireNamespace("tibble", quietly = TRUE)) {
#   install.packages("tibble")  # Install tibble
# }
# 
# if (!requireNamespace("ggplot2", quietly = TRUE)) {
#   install.packages("ggplot2")  # Install ggplot2
# }

# Load libraries
library('phyloseq')
library('tibble')
library('ggplot2')


###Files input
otu_mat <- read.delim('Otu_table_uterusvagina_t100_scaled_3_rare_sourcename.tsv')
head(otu_mat)

samples <- read.delim("Meta_table_uterusvagina_t100.txt")
head(samples)

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("X.Classification") 

samples <- samples %>% 
  tibble::column_to_rownames("source_name") 

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

## Make relative abundance if needed
uteruseq_relative <- transform_sample_counts(uteruseq, function(x) x / sum(x)) 

# Save the OTU table as a TSV file
#otu_as_df <- t(as.data.frame(otu_table(uteruseq_relative)))
#write.table(otu_as_df, file = "relative_otu_table.tsv", sep = "\t", 
#            row.names = TRUE, col.names = NA, quote = FALSE)


# Filter the data based on the variable if needed
#H <- subset_samples(uteruseq_relative, Diagnosis == "H")


### Top n taxa
## Option 1: Extract top n taxa of all samples 
#TopNOTUs <- names(sort(taxa_sums(uteruseq_relative), TRUE)[1:20])  #change [1:n] as needed
#topnuteruseq <- prune_taxa(TopNOTUs, uteruseq)
#topnuteruseq_relative <- transform_sample_counts(topnuteruseq, function(x) x / sum(x))


## Option 2: Extract top n taxa at a single sample lvl and then merged
# Define a function to prune top N taxa for a single sample
prune_top_taxa_per_sample <- function(physeq_obj, top_n = 10) {
  otu_matrix <- otu_table(physeq_obj)
  sample_names <- colnames(otu_matrix)  # Extract sample names
  
  # Create a list to store pruned objects for each sample
  pruned_list <- lapply(sample_names, function(sample) {
    # Extract data for the specific sample
    sample_data <- otu_matrix[, sample, drop = FALSE]
    
    # Identify top N taxa for the sample
    top_taxa <- rownames(sort(sample_data, decreasing = TRUE)[1:top_n])
    
    # Prune the phyloseq object to keep only these top taxa
    pruned_sample <- prune_taxa(top_taxa, prune_samples(sample, physeq_obj))
    
    return(pruned_sample)
  })
  
  names(pruned_list) <- sample_names  # Name each object in the list by the sample
  return(pruned_list)
}

# Apply the function to your phyloseq object
# Combine all pruned phyloseq objects into one using do.call
# Transform to relative abundances

top_20_pruned_list <- prune_top_taxa_per_sample(physeq_obj = uteruseq_relative, top_n = 20)
combined_20_uteruseq <- do.call(merge_phyloseq, top_20_pruned_list)
combined_20_uteruseq_relative <- transform_sample_counts(combined_20_uteruseq, function(x) x / sum(x))

top_10_pruned_list <- prune_top_taxa_per_sample(physeq_obj = uteruseq_relative, top_n = 10)
combined_10_uteruseq <- do.call(merge_phyloseq, top_10_pruned_list)
combined_10_uteruseq_relative <- transform_sample_counts(combined_10_uteruseq, function(x) x / sum(x))


### Set unique divergent color palette
palette65 <- c("#61a4d4",
               "#591f26",
               "#8789cd",
               "#6a752b",
               "#dda6e1",
               "#357140",
               "#e57ba4",
               "#1d2c19",
               "#e37756",
               "#132634",
               "#da9f7b",
               "#2f465f",
               "#a78f52",
               "#aa68a3",
               "#b0bdea",
               "#a75623",
               "#479193",
               "#99345c",
               "#325a58",
               "#d77a7d",
               "#587893",
               "#957fa1",
               "#785326",
               "#4d3925",
               "#d8afb9",
               "#674a63",
               "#798069",
               "#986866",
               "#da3b4e",
               "#a7d996",
               "#dd41a0",
               "#5325a6",
               "#bdcf38",
               "#a357df",
               "#65cf76",
               "#dd4425",
               "#e0bf38",
               "#51419d",
               "#85a42b",
               "#9c3196",
               "#55d5a6",
               "#5d2edf",
               "#458631",
               "#abc768",
               "#351964",
               "#e0892e",
               "#87af74",
               "#5a64e2",
               "#57cfdf",
               "#d93b72",
               "#9d75d6",
               "#4ea588",
               "#a88926",
               "#44538f",
               "#d9be6b",
               "#281e48",
               "#dd73d4",
               "#9fd4ab",
               "#753274",
               "#5082de",
               "#521845",
               "#c4c09a",
               "#321a27",
               "#9dc3ca",
               "#8b2c24")


palette38  <- c("#a75623", "#458631", "#c4c09a", "#d77a7d", "#325a58",
                "#2f465f", "#57cfdf", "#957fa1", "#674a63", "#591f26",
                "#986866", "#da3b4e", "#87af74", "#9c3196", "#d8afb9",
                "#1d2c19", "#e37756", "#e57ba4", "#b0bdea", "#da9f7b",
                "#4ea588", "#357140", "#d9be6b", "#dd41a0", "#51419d",
                
                "#9d75d6", "#521845", "#aa68a3", "#5a64e2", "#dda6e1",
                "#5325a6", "#bdcf38", "#e0bf38", "#9fd4ab", "#9dc3ca",
                "#753274", "#55d5a6", "#99345c")


paletteX <- sample(palette65) #Randomization of colors
taxanames <- rownames(otu_table(combined_10_uteruseq_relative))
unique_colors <- rep(paletteX, length.out = 80) #change length.out value for current number of taxa, colors are looped
unique_colors_df <- data.frame(Taxa = taxanames, Color = unique_colors)

# Adjust a palette to different taxa ranks
# Extract the part after 's__' using sub()
unique_colors_df$Species <- sub(".*s__", "", unique_colors_df$Taxa)
unique_colors_df$Genus <- sub(".*g__(.*?)\\|s__.*", "\\1", unique_colors_df$Taxa)
unique_colors_df_vector_species <- setNames(unique_colors_df$Color, unique_colors_df$Species)
unique_colors_df_vector_genus <- setNames(unique_colors_df$Color, unique_colors_df$Genus)

# Save colors df if needed
write.table(unique_colors_df, file = "unique_colors_df.tsv", sep = "\t", 
            row.names = TRUE, col.names = NA, quote = FALSE)


### Set sample order for plotting if needed
sample_names <- rownames(sample_data(combined_10_uteruseq_relative))
ordered_samples1 <- c("uterus_546", "uterus_629", "uterus_636", "uterus_633", 
                     "uterus_603", "uterus_634", "uterus_621", "uterus_700", 
                     "uterus_281", "uterus_270", "uterus_654", "uterus_709",
                     "uterus_594", "uterus_595", "uterus_597", "uterus_599",
                     "uterus_793", "uterus_799", "uterus_542",  "uterus_632",
                     "uterus_630", "uterus_966")

### Plot top n taxa
plot_bar(combined_10_uteruseq_relative, y = "Abundance", fill = "Genus") +
  facet_grid(~Microbiome_profile, scale="free") +
  xlab("Sample name") +
  ylab("Relative abundance") +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_fill_manual(values=unique_colors_df_vector_genus) 
  scale_fill_manual(values=palette38)
  #theme(legend.position = "none") +
  #scale_x_discrete(limits = ordered_samples1)

plot_bar(combined_10_uteruseq_relative, y = "Abundance", fill = "Species") +
  facet_grid(~Microbiome_profile, scale="free") +
  xlab("Sample name") +
  ylab("Relative abundance") +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=unique_colors_df_vector_species) 
  #theme(legend.position = "none") 
 # scale_x_discrete(limits = ordered_samples1)

plot_bar(combined_10_uteruseq_relative, y = "Abundance", fill = "Genus") +
  facet_grid(~Menopause_status~BMI, scale="free") +
  xlab("Sample name") +
  ylab("Relative abundance") +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=palette38) +
  theme(legend.position = "none")

# Save plots manually
# Save the OTU table as a TSV file
otu_as_df <- t(as.data.frame(otu_table(combined_10_uteruseq_relative)))
write.table(otu_as_df, file = "uterus_top10_uteruseq_relative_otu_table.tsv", sep = "\t", 
            row.names = TRUE, col.names = NA, quote = FALSE)

otu_as_df <- t(as.data.frame(otu_table(uteruseq_relative)))
write.table(otu_as_df, file = "uteruseq_relative_otu_table.tsv", sep = "\t", 
            row.names = TRUE, col.names = NA, quote = FALSE)



### UTERUS AND VAGINA PLOT_______________________________________________________
## Set more colors

palette42  <- c("#a75623", "#458631", "#c4c09a", "#d77a7d", "#325a58",
                "#2f465f", "#57cfdf", "#957fa1", "#674a63", "#591f26",
                "#986866", "#da3b4e", "#87af74", "#9c3196", "#d8afb9",
                "#1d2c19", "#e37756", "#321a27", "#e57ba4", "#b0bdea", "#d9be6b",
                "#4ea588", "#357140", "#d9be6b", "#dd41a0", "#51419d",
                
                "#9d75d6", "#521845", "#aa68a3", "#5a64e2", "#dda6e1",
                "#5325a6", "#bdcf38", "#e0bf38", "#9fd4ab", "#2f465f", "#9dc3ca",
                "#753274", "#55d5a6", "#99345c", "#9dc3ca", "#8b2c24")


"#da9f7b"
paletteX <- sample(palette65) #Randomization of colors
taxanames <- rownames(otu_table(combined_10_uteruseq_relative))
unique_colors <- rep(paletteX, length.out = 98) #change length.out value for current number of taxa, colors are looped
unique_colors_df <- data.frame(Taxa = taxanames, Color = unique_colors)

# Adjust a palette to different taxa ranks
# Extract the part after 's__' using sub()
unique_colors_df$Species <- sub(".*s__", "", unique_colors_df$Taxa)
unique_colors_df$Genus <- sub(".*g__(.*?)\\|s__.*", "\\1", unique_colors_df$Taxa)
unique_colors_df_vector_species <- setNames(unique_colors_df$Color, unique_colors_df$Species)
unique_colors_df_vector_genus <- setNames(unique_colors_df$Color, unique_colors_df$Genus)

# Save colors df if needed
write.table(unique_colors_df, file = "unique_colors_df.tsv", sep = "\t", 
            row.names = TRUE, col.names = NA, quote = FALSE)

## Set order manually if needed
ordered_samples2 <- c("patient_213","patient_220",
                      "patient_217","patient_190",
                      "patient_218","patient_210",
                      "patient_236","patient_036",
                      "patient_051","patient_207",
                      "patient_239","patient_093",
                      "patient_184","patient_269",
                      "patient_274","patient_216",
                      "patient_214","patient_384")

ordered_samples3 <- c("uterus_629","vagina_694",
                      "uterus_636","vagina_681", 
                      "uterus_633","vagina_652", 
                      "uterus_603","vagina_642",
                      "uterus_634","vagina_692",
                      "uterus_621", "vagina_622",
                      "uterus_700", "vagina_701",
                      "uterus_281", "vagina_382",
                      "uterus_270", "vagina_358",
                      "uterus_654", "vagina_616",
                      "uterus_709","vagina_710",
                      "uterus_594","vagina_645",
                      "uterus_597", "vagina_651",
                      "uterus_793", "vagina_794",
                      "uterus_799","vagina_800",
                      "uterus_632","vagina_618",
                      "uterus_630","vagina_696",
                      "uterus_966","vagina_967")

## Plot top n taxa
plot_bar(combined_10_uteruseq_relative, y = "Abundance", fill = "Genus") +
  facet_grid(~Patient_code, scale="free") +
  xlab("Sample name") +
  ylab("Relative abundance") +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=palette42) 
  #theme(legend.position = "none") 
  #scale_x_discrete(limits = ordered_samples3)

plot_bar(combined_10_uteruseq_relative, y = "Abundance", fill = "Species") +
  facet_grid(~Patient_code, scale="free") +
  xlab("Sample name") +
  ylab("Relative abundance") +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=unique_colors_df_vector_species) +
  theme(legend.position = "none") 
  #scale_x_discrete(limits = ordered_samples3)


# Save plots manually
# Save the OTU table as a TSV file
otu_as_df <- t(as.data.frame(otu_table(combined_10_uteruseq_relative)))
write.table(otu_as_df, file = "uterusvagina_top10_uteruseq_relative_otu_table.tsv", sep = "\t", 
            row.names = TRUE, col.names = NA, quote = FALSE)

otu_as_df <- t(as.data.frame(otu_table(uteruseq_relative)))
write.table(otu_as_df, file = "uterusvagina_uteruseq_relative_otu_table.tsv", sep = "\t", 
            row.names = TRUE, col.names = NA, quote = FALSE)


##Plot all microbiome profiles separately
# Filter the data based on the variable if needed
profile1 <- subset_samples(combined_10_uteruseq_relative, Microbiome_profile == "1")
profile2 <- subset_samples(combined_10_uteruseq_relative, Microbiome_profile == "2")
profile3 <- subset_samples(combined_10_uteruseq_relative, Microbiome_profile == "3")
profile4 <- subset_samples(combined_10_uteruseq_relative, Microbiome_profile == "4")


## Plot top n taxa
plot_bar(profile4, y = "Abundance", fill = "Genus") +
  facet_grid(~Microbiome_profile~Patient_code, scale="free") +
  xlab("Sample name") +
  ylab("Relative abundance") +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=palette42) +
  theme(legend.position = "none") 
  #scale_x_discrete(limits = ordered_samples3)

                                                         
###MERGED PROFILES ON ONE CHART________________________________
###Files input
otu_mat <- read.delim('Otu_table_uterus_merged_profiles.tsv')
head(otu_mat)

samples <- read.delim("Meta_table_uterus_merged_profiles.txt")
head(samples)

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("X.Classification") 

samples <- samples %>% 
  tibble::column_to_rownames("X.NAME") 

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

## Make relative abundance if needed
uteruseq_relative <- transform_sample_counts(uteruseq, function(x) x / sum(x)) 

top_5_pruned_list <- prune_top_taxa_per_sample(physeq_obj = uteruseq_relative, top_n = 5)
combined_5_uteruseq <- do.call(merge_phyloseq, top_5_pruned_list)
combined_5_uteruseq_relative <- transform_sample_counts(combined_5_uteruseq, function(x) x / sum(x))

## Set colors
palette17a  <- c("#87af74", "#d77a7d", "#957fa1", "#674a63", "#da3b4e",
                 "#b0bdea", "#d77a7d", "#51419d", "#55d5a6", "#51419d",
                 "#51419d", "#51419d", "#51419d", "#9d75d6", "#5325a6",
                 "#9fd4ab", "#55d5a6")

palette17b  <- c("#87af74", "#d77a7d", "#957fa1", "#674a63", "#da3b4e",
                 "#b0bdea", "#99345c", "#51419d", "#55d5a6", "#9c3196",
                 "#d9be6b", "#5a64e2", "#e57ba4", "#9d75d6", "#da9f7b",
                 "#325a58", "#1d2c19")

## Plot results
plot_bar(combined_5_uteruseq_relative, y = "Abundance", fill = "Species") +
  facet_grid(~Microbiome_profile, scale="free") +
  xlab("Microbiome_profile") +
  ylab("Relative abundance") +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  #scale_fill_manual(values=unique_colors_df_vector_genus) 
  scale_fill_manual(values=palette17)
#theme(legend.position = "none") +
#scale_x_discrete(limits = ordered_samples1)

plot_bar(combined_5_uteruseq_relative, y = "Abundance", fill = "Species") +
  facet_grid(~Microbiome_profile, scale="free") +
  xlab("Microbiome_profile") +
  ylab("Relative abundance") +
  geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=palette17b) 
#theme(legend.position = "none") 
# scale_x_discrete(limits = ordered_samples1)
