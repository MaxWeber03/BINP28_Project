R.version.string 
# R version 4.5.2 (2025-10-31)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("SNPRelate")
library(SNPRelate)
sessionInfo()
# SNPRelate_1.44.0 gdsfmt_1.46.0


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

# run pruning
spn_pruned = snpgdsLDpruning(snp_data, ld.threshold=0.2)
# why 0.2? The manual recommends to try different values. What am I looking for?

# SNP pruning based on LD:
#   Excluding 50,188 SNPs on non-autosomes
# Excluding 710 SNPs (monomorphic: TRUE, MAF: 0.005, missing rate: 0.01)
# # of samples: 15
# # of SNPs: 61,431
# using 1 thread/core
# sliding window: 500,000 basepairs, Inf SNPs
# |LD| threshold: 0.2
# method: composite
# Chrom 5: |====================|====================|
#   1.12%, 696 / 62,141 (Sun Feb 15 18:08:05 2026)
# 696 markers are selected in total.

str(spn_pruned)
# only chr5 remains.


# PCA ---------------------------------------------------------------------
snpset.id <- unlist(unname(spn_pruned))
pca <- snpgdsPCA(snp_data, snp.id=snpset.id, num.thread=2)


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

# create table with eigenvectors
eigenvec_table1234 = data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
print(eigenvec_table1234)

# Draw 1 and 2
plot(x = eigenvec_table1234$EV2, 
     y = eigenvec_table1234$EV1, 
     col=as.integer(eigenvec_table1234$pop), 
     xlab="eigenvector 2", 
     ylab="eigenvector 1")
legend("bottomleft", 
       legend=levels(eigenvec_table1234$pop), 
       pch="o", col=1:nlevels(eigenvec_table1234$pop))


# Draw 3 and 4
plot(x = eigenvec_table1234$EV4, 
     y = eigenvec_table1234$EV3, 
     col=as.integer(eigenvec_table1234$pop), 
     xlab="eigenvector 4", 
     ylab="eigenvector 3")
legend("bottom", 
       legend=levels(eigenvec_table1234$pop), 
       pch="o", col=1:nlevels(eigenvec_table1234$pop))


# Calculate explained variance per eigenvector ----------------------------
# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
# 11.20  8.80  8.57  7.96  7.42  7.38



# Correlation -------------------------------------------------------------


