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
