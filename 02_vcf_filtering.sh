###########
# Step 0: Install bcftools via conda
###########

# Install vcftools from bioconda (current version 0.1.17)
conda create -n bcftools bcftools=1.23
conda install matplotlib
# necessary for plotting with bcftools


# Activate env
conda activate bcftools

# Confirm version
bcftools --version
# 1.23 using htslib 1.23
conda list matplotlib
# 3.9.2

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

grep -v ^\## 02_vcf_filtered/outgroup_removed.vcf | less -S

# The QUAL values are not changed after removing Naxos2 and hence not valid anymore.
# There will be variants now, that do not really exist, since they were just introduced by the outgroup. This can be filtered with AC.

# For the AN we should habe AN=30 if all samples have been genotyped. This is a very strong filter, but this would ensure, that we actually have data on all of our samples for the site.

# Before filtering on AC/AN we can filter on each sample on the genotype quality and/or the number of reads.

# Masking samples for some sites, will reduce the AN/AC. We can then filter based on AN/AC in a following step.

###########
# 2.1 Filtering per sample sequencing depth
###########

# We want to filter out very numbers of reads, as these are likely artifacts, and very low numbers of reads, since these have not so much confidence. To have an idea what few or many is in our data, we should look at the distribution of sequenceing depth.

# We can get summary stats from bcftools:
mkdir 02_vcf_filtered/stats_plots

bcftools stats -s - 02_vcf_filtered/outgroup_removed.vcf > 02_vcf_filtered/stats_plots/stats.vchk
# this generates stats on the vcf file -s - includes all the samples

# create stats
plot-vcfstats -p 02_vcf_filtered/stats_plots/ 02_vcf_filtered/stats_plots/stats.vchk

# Let's filter an > 30 (maximum one sample not genotyped)
# and AC > 3 (about 10% of the AN must be different from the reference)
# todo: filter for sequencing depth
# bcftools view -i 'INFO/AN>30 && INFO/AC>3' \



