# Max Weber BINP28 Project
# Try to get LD decay plot

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


# Open file ---------------------------------------------------------------

snp_data = snpgdsOpen("./03_gds_pruning/converted.gds")
str(snp_data)

# LD decay plot -----------------------------------------------------------
# Create LD decay plot to find a sensible ld.threshold 
# SNPRelate offers different methods for calculating the LD between two Variants, and for LD pruning
?snpgdsLDMat
?snpgdsLDpruning

# method=c("composite", "r", "dprime", "corr", "cov")

# Take sample of variants (using all them would be too big
ld_decay_snp.id = sample(snpgdsSNPList(snp_data)$snp.id, 1000)
length(ld_decay_snp.id) #1000 snp.ids

# next we can use snpgdsLDMat to get the LD for pairs of variants
ld_result = snpgdsLDMat(
  gdsobj = snp_data,
  snp.id = ld_decay_snp.id,
  method = "r") 
# selects our data (snp_data), the sampled snp.ids and returns the LD as r

ld_r = ld_result$LD
hist(ld_r)

# the LD is almost normally distributed, which is not what I expect.
# By looking at LD decay plots online, I would expect more data points with lower ld values

# get distances (bp) between the variants
snpinfo <- snpgdsSNPList(snp_data)
snp.pos <- snpinfo$pos[ld_decay_snp.id]
distances <- as.vector(abs(outer(snp.pos, snp.pos, "-")))
# outer() takes the positions of the snp and makes a matrix of all pairwise distances
# abs() then gets the absolute value of each distance
length(distances)

# if the distance is zero we are comparing a variant against itself
keep_trues = distances != 0

# exclude these through the false/true values
distances = distances[keep_trues]
ld_r = ld_r[keep_trues]
length(distances)
str(distances)
length(ld_r)
str(ld_r)
# we now have removed some of the values, and the two vectors have the same length

# make df
ld_decay_df = data.frame(
  distance = distances,
  ld_r = ld_r
)

# sample from it to have a plot with not too many points.
ld_decay_df_sampled = slice_sample(
  .data = ld_decay_df,
  n = 1000
)
nrow(ld_decay_df_sampled)

plot(x = ld_decay_df_sampled$distance, y = ld_decay_df_sampled$ld_r)

# In this plot, there is not clear trend. The number of samples from the variants 
# and from the the ld_decay_df had been altered, but the pattern of the plot had not changed.
# I also checked the methods "corr" "dprime" and "composite". The results were different, but none resulted in a 
# plot that matched the expected pattern.
# The plots of different methods can be made by replacing the method in the function call snpgdsLDMat() above.

# LD Pruning -------------------------------------------------------------

# the pruning uses some randomeness, for reproducibility, we can set a seed
set.seed(seed = 1)

spn_pruned = snpgdsLDpruning(snp_data, ld.threshold=0.2, autosome.only = FALSE, method="r")
str(spn_pruned)

# By looking at the plot, I would expect more then 2% remaning. However, the threshold of 0.2 was also checked with plink2
# with one of the TAs and yielded 9000 markers there. That is more then the 2000 markers left here, but that could be due to 
# differences in the window sizes, and still points towards a high fraction of my markers (>90%) not passing the threshold 0.2 in
# both SNPRelate and Plink.