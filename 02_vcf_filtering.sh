# Install vcftools from bioconda (current version 0.1.17)
conda create -n vcftools -c bioconda vcftools=0.1.17

# Activate env
conda activate vcftools

# Confirm version
vcftools --version
# VCFtools (0.1.17)

man vcftools
