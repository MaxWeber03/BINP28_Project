R.version.string 
# R version 4.5.2 (2025-10-31)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPRelate")
library(SNPRelate)
