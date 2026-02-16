# Max Weber BINP28 Project
# Import filtered vcf and analyse with SNPRelate to get PCA

# Packages and Version Info -----------------------------------------------

R.version.string 
# R version 4.5.2 (2025-10-31)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("SNPRelate")
library(SNPRelate)
library(tidyverse)
sessionInfo()
# SNPRelate_1.44.0 gdsfmt_1.46.0 tidyverse_2.0.0
 

# Open GDS file -----------------------------------------------------------

snp_data = snpgdsOpen("./03_gds_pruning/converted.gds")
str(snp_data)

# LD decay plot -----------------------------------------------------------
# Create LD decay plot to find a sensible ld.threshold (r²)

# Take sample of Variants (using all them would be too big
ld_decay_snp.id = sample(snpgdsSNPList(snp_data)$snp.id, 500)
length(ld_decay_snp.id) #1000 snp.ids

# next we can use snpgdsLDMat to get the LD for pairs of variants
ld_result = snpgdsLDMat(gdsobj = snp_data, snp.id = ld_decay_snp.id, method = "corr") 
# selects our data (snp_data), the sampled snp.ids and returns the LD as r (orther metrics are avilable too, 
# but r² is the one used downstream by snpgdsLDpruning() )

ld_composite = ld_result$LD
hist(ld_composite)


# get distances between the variants
snpinfo <- snpgdsSNPList(snp_data)
snp.pos <- snpinfo$pos[ld_decay_snp.id]
distances <- as.vector(abs(outer(snp.pos, snp.pos, "-")))
  # outer() takes the positions of the snp and makes a matrix of all pairwise distances
  # abs() then gets the absolute value of each distance
distances # distances
length(distances)

# if the distance is zero we are comparing a variant against itself
keep_trues = distances != 0

# exclude these through the false/true values
distances = distances[keep_trues]
ld_composite = ld_composite[keep_trues]
length(distances)
str(distances)
length(ld_composite)
str(ld_composite)
# we now have removed some of the values, and the two vectors have the same length

# make df
ld_decay_df = data.frame(
  distance = distances,
  ld_composite = ld_composite
)

ld_decay_df_sampled = slice_sample(
  .data = ld_decay_df,
  n = 1000
)
nrow(ld_decay_df_sampled)

plot(x = ld_decay_df_sampled$distance, y = ld_decay_df_sampled$ld_composite)

# LD Pruning -------------------------------------------------------------

# the pruning uses some randomeness, for reproducibility, we can set a seed
set.seed(seed = 1)

spn_pruned = snpgdsLDpruning(snp_data, ld.threshold=0.2, autosome.only = FALSE, method="r")
str(spn_pruned)

# PCA ---------------------------------------------------------------------

snpset.id <- unlist(unname(spn_pruned))
pca <- snpgdsPCA(snp_data, snp.id=snpset.id, algorithm = "exact")
# algorithm can be randomized (faster, causes small deviations) or exact. 
# Since our dataset it not that big, we will go with exact for better repdorcuibility


# PCA Plotting ------------------------------------------------------------
# SNPRelate would like to have population info, let's make a vector for that
sample.id = read.gdsn(index.gdsn(snp_data, "sample.id"))
# make a new vector of the populations (I will name according to their naming, since we do not know the populations)
pop_code = c("8N", "8N","8N","8N","8N","K","K","K","K","K", "Lesina", "Lesina", "Lesina", "Lesina", "Lesina")
# confirm the order matches
print(
  cbind(
    sample.id, 
    pop_code)) 

# create table with eigenvectors 1-4, sample names and populations
eigenvec_table1234 = data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4])
print(eigenvec_table1234)

#SNPRelate calls the axis eigenvector 1 and 2, but it is the same as pc 1 and 2
# https://support.bioconductor.org/p/119389/

# Draw 1 and 2
# plot(x = eigenvec_table1234$EV2, 
#      y = eigenvec_table1234$EV1, 
#      col=as.integer(eigenvec_table1234$pop), 
#      xlab="PC 2", 
#      ylab="PC 1")
# legend("bottomleft", 
#        legend=levels(eigenvec_table1234$pop), 
#        pch="o", col=1:nlevels(eigenvec_table1234$pop))
# 
# 
# # Draw 3 and 4
# plot(x = eigenvec_table1234$EV4, 
#      y = eigenvec_table1234$EV3, 
#      col=as.integer(eigenvec_table1234$pop), 
#      xlab="PC 4", 
#      ylab="PC 3")
# legend("bottom", 
#        legend=levels(eigenvec_table1234$pop), 
#        pch="o", col=1:nlevels(eigenvec_table1234$pop))


# Rebuild the graph in ggplot (nicer, easier export) ----------------------
# PC1-2
pc12 = ggplot(data = eigenvec_table1234,
       aes(x = EV1, y = EV2, color = pop)) +
  geom_point() +
  labs(
    x = "PC1",
    y = "PC2",
    color = "Population"
  ) +
  theme(legend.position = "bottom")

pc12

# PC3-4
pc34 = ggplot(data = eigenvec_table1234,
       aes(x = EV3, y = EV4, color = pop)) +
  geom_point() +
  labs(
    x = "PC3",
    y = "PC4",
    color = "Population"
  ) +
  theme(legend.position = "bottom")
  
pc34

ggsave(filename = "./04_pca/pc12.png", plot = pc12, dpi = 300, width = 5, height = 4)
ggsave(filename = "./04_pca/pc34.png", plot = pc34, dpi = 300, width = 5, height = 4)


# Calculate explained variance per eigenvector ----------------------------

# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
# 11.19  8.67  8.62  7.83  7.42  7.31


# Make Elbow/Scree Plot ---------------------------------------------------

elbow_plot_data = 
  data.frame(
    variance = pc.percent,
    PC = seq(1:length(pc.percent))
  )

ggplot(data = elbow_plot_data, aes(x = PC, y = variance)) +
  labs(x = "PC", y = "Variance accounted for per PC") +
  geom_point()

# SNP Loadings - How much does each SNP contribute to a PC ----------------

snp_loadings = snpgdsPCASNPLoading(
  pcaobj = pca,
  gdsobj = snp_data
)

names(snp_loadings)
dim(snp_loadings$snploading) # SNP loadings, or SNP eigenvectors
# View(snp_loadings$snploading)
# we now have a df of loadings with PC in rows and snp in columns and the loading from the PCA as values

plot(snp_loadings$snploading[1,], type="h", ylab="PC 1")

# How can this be connected back to the SNPs? What value/unit does the load have?
length(snp_loadings$snp.id)
# We have a SNP ID for each loading

# So we should be able to find that number in our an_ac_filtered file,if we have IDs in that file.
# As of now, our file has no IDs (just ".") so the ID we have here is internal from SNPRelate.
# It should be possible make this work by giving proper IDs to the vcf file (maybe with bcftools?)
# Alternativly snpgdsSNPList() gives out the snp.id chromosome etc for each snp
snpgdsSNPList(snp_data)
