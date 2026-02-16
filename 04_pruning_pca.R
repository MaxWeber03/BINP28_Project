# Max Weber BINP28 Project
# Import gdp and analyse with SNPRelate to get PCA

# Packages and Version Info -----------------------------------------------

R.version.string 
# R version 4.5.2 (2025-10-31)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("SNPRelate")
# install.packages("tidyverse)
library(SNPRelate)
library(tidyverse)
sessionInfo()
# SNPRelate_1.44.0 gdsfmt_1.46.0 tidyverse_2.0.0
 

# Open GDS file -----------------------------------------------------------

snp_data = snpgdsOpen("./03_gds_pruning/converted.gds")
str(snp_data)

# LD Pruning -------------------------------------------------------------
# The LD decay plots (see 05_ld_decay.R) did not give any proper results.
# Hence I am now just sticking to the default value from the SNPRelate manual (0.2)

# the pruning uses some randomeness, for reproducibility, we can set a seed
set.seed(seed = 1)

spn_pruned = snpgdsLDpruning(snp_data, ld.threshold=0.2, autosome.only = FALSE, method="r")
# SNP pruning based on LD:
#   Excluding 1,067 SNPs (monomorphic: TRUE, MAF: 0.005, missing rate: 0.01)
# # of samples: 15
# # of SNPs: 111,262
# using 1 thread/core
# sliding window: 500,000 basepairs, Inf SNPs
# |LD| threshold: 0.2
# method: R
# Chrom 5: |====================|====================|
#   2.15%, 1,334 / 62,141 (Mon Feb 16 15:48:26 2026)
# Chrom Z: |====================|====================|
#   2.03%, 1,020 / 50,188 (Mon Feb 16 15:48:26 2026)
# 2,354 markers are selected in total.

# str(spn_pruned)

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

# Calculate explained variance per eigenvector ----------------------------

# variance proportion (%)
pc.percent <- pca$varprop*100
pc.percent
# sharp drop at PC 15


# Rebuild the graph in ggplot (nicer, easier export) ----------------------
# PC1-2
pc12 = ggplot(data = eigenvec_table1234,
       aes(x = EV1, y = EV2, color = pop)) +
  geom_point() +
  labs(
    x = paste0("PC1 (", round(pc.percent[1], 2), "% Variance)"),
    y = paste0("PC2 (", round(pc.percent[2], 2), "% Variance)"),
    color = "Population"
  ) +
  theme(legend.position = "bottom")

pc12

# PC3-4
pc34 = ggplot(data = eigenvec_table1234,
       aes(x = EV3, y = EV4, color = pop)) +
  geom_point() +
  labs(
    x = paste0("PC3 (", round(pc.percent[3], 2), "% Variance)"),
    y = paste0("PC4 (", round(pc.percent[4], 2), "% Variance)"),
    color = "Population"
  ) +
  theme(legend.position = "bottom")
  
pc34

ggsave(filename = "./04_pca/pc12.png", plot = pc12, dpi = 300, width = 5, height = 4)
ggsave(filename = "./04_pca/pc34.png", plot = pc34, dpi = 300, width = 5, height = 4)


# Make Elbow/Scree Plot ---------------------------------------------------

elbow_plot_data = 
  data.frame(
    variance = pc.percent,
    PC = seq(1:length(pc.percent))
  )

elbow_plot = ggplot(data = elbow_plot_data, aes(x = PC, y = variance)) +
  labs(x = "PC", y = "Variance accounted for per PC (%)") +
  geom_point()

elbow_plot

ggsave(filename = "./04_pca/elbow_plot.png", plot = elbow_plot, dpi = 300, width = 5, height = 4)

# First 4 components explain variance, but the amout of variance explained after the first 4 does not drop of clearly, but very slowly
# Not the typical elbow shape

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
