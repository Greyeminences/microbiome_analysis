# Install the required packages if they are not already installed
packages <- c("microbiome", "DirichletMultinomial", "reshape2", "magrittr", "dplyr")

# Install missing packages
#install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)

# Load the libraries
lapply(packages, library, character.only = TRUE)
library(DirichletMultinomial)

niche <- 'uterus'
output_folder <- '~/Huterus_again/again'

###Files input
otu_mat <- read.delim('Otu_table_all_t0_none_1_rare.tsv')
head(otu_mat)

samples <- read.delim("Meta_table_all_t0.txt")
head(samples)

otu_mat <- otu_mat %>%
  group_by(Classification) %>%  
  summarise(across(everything(), sum, na.rm = TRUE))

otu_mat <- otu_mat %>%
  tibble::column_to_rownames("Classification") 

samples <- samples %>% 
  tibble::column_to_rownames("sample_name") 

otu_mat <- as.matrix(otu_mat)

# Extract the rownames (GTDB long names) from the OTU table
gtdb_long_names <- rownames(otu_mat)

# Split the GTDB long names into separate taxonomic ranks by the ';' delimiter
split_taxa <- strsplit(gtdb_long_names, ";") %>%
  lapply(function(x) sub("^[a-z]__", "", x))

# Find the maximum length of any taxonomy vector
max_length <- max(sapply(split_taxa, length))

# Pad all vectors to this length with NA values
split_taxa <- lapply(split_taxa, function(x) {
  length(x) <- max_length
  x
})

# Now you can safely rbind
taxonomy <- do.call(rbind, split_taxa) %>%
  as.data.frame()

# Assign column names to represent taxonomic ranks
colnames(taxonomy) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")[1:ncol(taxonomy)]

# Set rownames (corresponding to OTUs) to match the OTU table
rownames(taxonomy) <- rownames(otu_mat)

# Convert to matrix for phyloseq
tax_matrix <- as.matrix(taxonomy)

# Check the first few rows of the generated taxonomy table
head(tax_matrix)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
samples = sample_data(samples)

cancerseq <- phyloseq(OTU, TAX, samples)
cancerseq

###change names
# Step 1: Extract sample data and ensure correct order
sd <- as.data.frame(sample_data(cancerseq))
sd <- sd[sample_names(cancerseq), ]

# Step 2: Get new names from source_name column
new_names <- sd$source_name

# Step 3: Safety checks
stopifnot(length(new_names) == nsamples(cancerseq))
stopifnot(!anyDuplicated(new_names))
stopifnot(all(rownames(sd) == sample_names(cancerseq)))

# Step 4: Assign new sample names
sample_names(cancerseq) <- new_names

###
# Filter the data based on the variable if needed
#H <- subset_samples(cancerseq_relative, Diagnosis == "H")
uterus <- subset_samples(cancerseq, Source == "uterus")
vagina <- subset_samples(cancerseq, Source == "vagina")

pseq <- uterus

# Identify core taxa using presence/absence (prevalence-based)
pseq.pa <- transform_sample_counts(pseq, function(x) ifelse(x > 0, 1, 0))
taxa <- core_members(pseq.pa, detection = 0, prevalence = 0.1) # 10% prevalence
pseq.core <- prune_taxa(taxa, pseq)

# Extract RAW COUNTS (not relative abundances)
count <- as.matrix(t(abundances(pseq.core)))
count <- ceiling(as.matrix(count))

# Verify integer values
print(head(count[, 1:5])) # Should show whole numbers

# Fit DMM
fit <- lapply(1:7, dmn, count = count, verbose = TRUE)


#Check model fit with different number of mixture components using standard information criteria
lplc <- base::sapply(fit, DirichletMultinomial::laplace) # AIC / BIC / Laplace
aic  <- base::sapply(fit, DirichletMultinomial::AIC) # AIC / BIC / Laplace
bic  <- base::sapply(fit, DirichletMultinomial::BIC) # AIC / BIC / Laplace

plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)
legend("topright", legend=c("Laplace", "AIC", "BIC"), lty=1:3)



#Pick the optimal model
best <- fit[[which.min(unlist(lplc))]]

#Mixture parameters pi and theta
mixturewt(best)


#Sample-component assignments
ass <- apply(mixture(best), 1, which.max)


#Contribution of each taxonomic group to each component
for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU)))
    # Only show the most important drivers
    #filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}


# Assuming 'best' is your selected DMN model object
fit_matrix <- fitted(best)  # taxa x components matrix of parameter estimates

# For each component k, extract and plot top taxa
for (k in seq_len(ncol(fit_matrix))) {
  df <- melt(fit_matrix[, k, drop=FALSE])
  colnames(df) <- c("OTU", "cluster", "value")
  
  # Arrange taxa by importance in component k
  df <- df %>%
    arrange(desc(value)) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU)))
  
  df <- df %>%
    mutate(Taxon = str_trim(sapply(str_split(OTU, ";"), function(x) tail(x, 1)))) %>%
    arrange(desc(value)) %>%
    mutate(Taxon = factor(Taxon, levels = unique(Taxon)))
  
  # Optionally filter top taxa, e.g., top 20
  top_df <- head(df, 20)
  
  # Plot
  p <- ggplot(top_df, aes(x = Taxon, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top taxa driving component", k))
  print(p)
}
