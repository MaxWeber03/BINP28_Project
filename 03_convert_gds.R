# Max Weber BINP28 Project
# Import filtered vcf and analyse with SNPRelate to get PCA

# Packages and Version Info -----------------------------------------------

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