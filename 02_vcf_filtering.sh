###########
# Step 0: Install bcftools via conda
###########

# Install vcftools from bioconda (current version 0.1.17)
conda create -n bcftools bcftools=1.23

# Activate env
conda activate bcftools

# Confirm version
bcftools --version
# 1.23 using htslib 1.23

# Manual
man bcftools

# Directory for filtered data
mkdir 02_vcf_filtered


###########
# Step 1: Remove Outgroup "Naxos2"
###########

# Removing a sample changes the true INFO values (AC, AN). These may be automatically adjusted when used with bcftools view.
# --sample ^<sample> will exclude a sample

bcftools view --samples ^Naxos2 01_raw/ProjTaxa.vcf \
    --output 02_vcf_filtered/outgroup_removed.vcf

# Check raw file
cat 01_raw/ProjTaxa.vcf | grep -v ^\## | less -S
# Check file without Naxos2
cat 02_vcf_filtered/outgroup_removed.vcf | grep -v ^\## | less -S
# AC and AN were actually adjusted

###########
# Step 2: Filtering
###########

