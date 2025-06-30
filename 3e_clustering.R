
# install.packages("factoextra")
library("cluster")
library("factoextra")
library("magrittr")
#install.packages("pheatmap")
library(pheatmap)


### Load otu file and create phyloseq object
unifrac_dist <- read.delim('otu_uterus_t0_1_rare_unifracc.txt', sep = '\t')

unifrac_dist <- unifrac_dist %>%
  tibble::column_to_rownames("X") 
unifrac_mat <- as.matrix(unifrac_dist)

#samples <- read.delim("Meta_table_uterus_t0.txt")
samples <- read.delim("metadata_u_4_clusters.txt")
head(samples)

samples <- samples %>% 
  tibble::column_to_rownames("sample_name") 

# Create a lookup vector
name_map <- samples$source_name
names(name_map) <- rownames(samples)

# Replace row and column names in the matrix
rownames(unifrac_mat) <- name_map[rownames(unifrac_mat)]
colnames(unifrac_mat) <- name_map[colnames(unifrac_mat)]

rownames(samples) <- samples$source_name


#_________________________________
#### Determine nb of clusters first

### Average Silhouette Width (ASW)
library(cluster)

# Try multiple k values
sil_widths <- sapply(2:10, function(k) {
  pam_fit <- pam(unifrac_mat, k = k)
  mean(silhouette(pam_fit)[, "sil_width"])
})

# Plot
plot(2:10, sil_widths, type = "b", pch = 19,
     xlab = "Number of clusters", ylab = "Average silhouette width")


### Gap statistic
library(cluster)
library(factoextra)

gap_stat <- clusGap(as.matrix(unifrac_mat), FUN = pam, K.max = 10, B = 100)
fviz_gap_stat(gap_stat)


# clValid – Internal & stability validation (similar to ConsensusClusterPlus)
library(clValid)

# Again, use Euclidean data (e.g., from PCoA)
coords <- cmdscale(unifrac_mat, k = 10)

# Evaluate clustering
cl_res <- clValid(coords, nClust = 2:6,
                  clMethods = c("hierarchical", "pam"),
                  validation = "internal")

summary(cl_res)
plot(cl_res)


#_________________________________
## k-medoids/PAM clustering
coords <- cmdscale(as.dist(unifrac_mat))
pam.res <- pam(coords, k = 3)

# Convert distance matrix to dist object
d <- as.dist(unifrac_mat)

# Compute silhouette info
sil_pam <- silhouette(pam.res$clustering, d)
sil_pam_df <- as.data.frame(sil_pam)

fviz_silhouette(pam.res)

sil_pam_df[sil_pam_df$sil_width < 0, ]

write.table(sil_pam_df, file = "silhouette_values_pam.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE)


## Compute hierarchical clustering
res.hc <- as.dist(unifrac_mat) %>%
  hclust(method = "ward.D2")  

fviz_dend(res.hc, k = 3, # Cut in four groups
          cex = 0.7, # label size
          k_colors = c( "#00AFBB", "#E7B800", "#753274"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)

hc_clusters <- cutree(res.hc, k = 3)
head(hc_clusters)

d <- as.dist(unifrac_mat)
hc <- hclust(d, method = "ward.D2") 
k <- 3
hc_clusters <- cutree(hc, k = k)
sil <- silhouette(hc_clusters, d)
fviz_silhouette(sil)
summary(sil)
sil_df <- as.data.frame(sil)
head(sil_df[order(sil_df$sil_width), ])



# Silhouette width of observations
sil <- res.hc$silinfo$widths[, 1:3]
sil
# Objects with negative silhouette
neg_sil_index <- which(sil[, 'sil_width'] < 0)
sil[neg_sil_index, , drop = FALSE]


# Or export to tab-delimited TXT
write.table(sil, file = "silhouette_values_hierarchical.txt", sep = "\t", quote = FALSE, row.names = TRUE)

#______________________
samples$Clusters_hierarchical <- cutree(res.hc, k = 3)
samples$Clusters_pam <- pam.res$clustering

#_________________________________
# Order the matrix by group
samples_ordered <- samples[order(samples$Clusters_hierarchical), ]
ordered_samples <- rownames(samples_ordered)

stopifnot(all(ordered_samples %in% rownames(unifrac_mat)))
stopifnot(all(ordered_samples == rownames(samples_ordered)))

unifrac_mat_ordered3 <- unifrac_mat[ordered_samples, ordered_samples]


annotation_df <- samples[, c('Clusters_hierarchical',
                             'Clusters_pam',
                             'Clusters_DMM')]

annotation_df <- samples[, c("Pregnancy_occurence",
                             "BMI",
                             "Menopause_status",
                             "Cycle_day",
                             "Age",
                             "Sample_collection_method")]

pheatmap(unifrac_mat_ordered3, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         display_numbers = FALSE,
         number_format = "%.4f",
         # main = "UniFrac distance matrix",
         annotation_row = annotation_df,
         annotation_col = annotation_df,
         #annotation_colors = ann_colors,
         border_color = NA,
         annotation_names_row = FALSE)

