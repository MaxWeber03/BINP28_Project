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
 

# Read .vcf and convert to GDS --------------------------------------------

snpgdsVCF2GDS(vcf.fn = "./02_vcf_filtered/an_ac_filtered.vcf",
              out.fn = "./03_gds_pruning/converted.gds",
              method = "copy.num.of.ref" # ncluding bi-allelic SNPs, multi-allelic SNPs, indels and structural variants
              )

snpgdsSummary(gds = "./03_gds_pruning/converted.gds")
# The file name: /home/max/OneDrive/Uni/Master_Lund/2._Semester/BINP28_Sequencing_Informatics_I/11_Project/03_gdp/converted.gdc 
# The total number of samples: 15 
# The total number of SNPs: 112329 
# SNP genotypes are stored in SNP-major mode (Sample X SNP).
# The number of valid samples: 15 
# The number of biallelic unique SNPs: 99971 


# LD Pruning --------------------------------------------------------------

# open data
snp_data = snpgdsOpen("./03_gds_pruning/converted.gds")
str(snp_data)

# the pruning uses some randomeness, for reproducibility, we can set a seed
set.seed(seed = 1)

# run pruning
spn_pruned = snpgdsLDpruning(snp_data, ld.threshold=0.2, autosome.only = FALSE)
# why 0.2? The manual recommends to try different values. What am I looking for?

str(spn_pruned)
# only chr5 remains when autosome.only = FALSE is set to TRUE (default)
# WIth FALSE we now have Chrom 5 ad Chrom Z


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

#SNPRelate calls the axises eigenvector 1 and 2, but it is the same as pc 1 and 2
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
