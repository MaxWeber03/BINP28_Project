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
# 2.1 Masking per sample sequencing depth
###########

# We want to filter out very numbers of reads, as these are likely artifacts, and very low numbers of reads, since these have not so much confidence. To have an idea what few or many is in our data, we should look at the distribution of sequenceing depth.

# We can get summary stats from bcftools:

bcftools stats -s - 02_vcf_filtered/outgroup_removed.vcf > 02_vcf_filtered/stats_outgroup_removed.vchk
# this generates stats on the vcf file -s - includes all the samples

# create stats
plot-vcfstats -v -p 02_vcf_filtered/stats_plots_outgroup_removed/ 02_vcf_filtered/stats_outgroup_removed.vchk
# -p is the output dir
# -v vectorized plots (clearer resolution)

# Interpretation of the Depth distribution plot:
# main peak at Depth of 7, declines around 10-15. Number of genotypes for low depth (<4) is small. The cummulitative number of genotypes increases strongly between 5 and 10-15.

# => From this, I am choosing to cut out genotypes with less then 5 and then 14 reads. This means we excluding around 10% of genotypes at the low end (likely to be errors) and about 10% at the high end (likely to be artifacts).

bcftools filter -S . -e 'FORMAT/DP<5 | FORMAT/DP>14' 02_vcf_filtered/outgroup_removed.vcf -o 02_vcf_filtered/dp_masked.vcf
# -S . let's us mask (replace with missing genotype) the values for samples for a variant, that did not pass the threshold
# -e applies the Threshold, FMT/DP means the FORMAT/DP field of each sample
# | => or, exclude if one of the two conditions is met

# confirm that AC and AN have been changed
grep -v ^\#\# 02_vcf_filtered/outgroup_removed.vcf | less -S
grep -v ^\#\# 02_vcf_filtered/dp_masked.vcf | less -S
# yes, some of the AN values have changed, so they appear as they were recalculated
# if there was more time, it would make sense to check this deeper

###########
# 2.1 Filtering AC and AN
###########

# Let's filter AN > 30 (all sites must be present in all samples, this could be done less harshely, but I lack time to mess around with it)
# and AC > 2 (excluse rare alternatives)

bcftools view -e 'INFO/AN<30 | INFO/AC<3' \
    02_vcf_filtered/dp_masked.vcf -o 02_vcf_filtered/an_ac_filtered.vcf


grep -v ^\#\# 02_vcf_filtered/an_ac_filtered.vcf | less -S
# filtering has worked

###########
# 3 Let's rerun the statistics plots to confirm that they have changed as expected
###########

bcftools stats -s - 02_vcf_filtered/an_ac_filtered.vcf > 02_vcf_filtered/after_filtering.vchk
# this generates stats on the vcf file -s - includes all the samples

# create stats
plot-vcfstats -v -p 02_vcf_filtered/stats_plots_after_filtering/ 02_vcf_filtered/after_filtering.vchk
# -p is the output dir
# -v vectorized plots (clearer resolution)

# ts/tv improved from 1.68 to 1.79, no more singletons (AC=1)
# the read depth plot has changed according to our cutoffs


# Delete data not needed anymore (replication should be no issue), to save disk space
rm -r 02_vcf_filtered/stats_plots_after_filtering
rm -r 02_vcf_filtered/stats_plots_outgroup_removed
rm 02_vcf_filtered/dp_masked.vcf
rm 02_vcf_filtered/*.vchk
rm 02_vcf_filtered/outgroup_removed.vcf

