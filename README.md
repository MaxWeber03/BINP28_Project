Max Weber BINP28 M. Sc. Bioinformatics Programme Lund University, Feb 2026

A Git Repository for this project can be found here: https://github.com/MaxWeber03/BINP28_Project.

# Instructions Summary
We have a .vcf file without much futher detail on the data or the species. I am tasked with "Genetic relationship/Multivariate relationship from a PCA". We do not have much information. 3.5 days are available including writing of a two pages report and a two slide/3 min presentation.

# Workflow Overview

The analysis here is divided into 4 major steps with respective folders and scripts. Running the scripts requires R and Conda to be installed. The installation commands for the used R packages is included in the R scripts as comments (remove # to run).

Software Versions:
 - conda 25.11.1
 - R version 4.5.2 (2025-10-31)
 - bcftools 1.23 with matplotlib 3.9.2
 - tidyverse R package 2.0.0
 - SNPRelate R package 1.44.0

0. __Report__: Not part of the analysis workflow, the directory 00_report contains the latex code to generate the report that needs to be submitted for this assignment. I know that having both the report (latex code) and the project code in the same repo is not a clean solution, but I want to avoid bloating my github with multiple repos for a small project.
1. __Download and unpacking__ of raw data into 01_raw directory, done by 01_data_download_extraction.sh
2. __Filtering of variants__: Done by 02_vcf_filtering.sh, output saved in 02_vcf_filtered directory. Removes outgrup sample from data. Filtering is done with bcftools by masking sites with less then 5 or more then 14 reads. These thresholds were determined from the histrogram of the read distribution produced by bcftool's plot-vcfstats function. The masked data was then filtered to not have AN below 30 (only variants without missing data, 15 diploid samples), and not less then 3 times the alterantive allele (AC > 2). No filtering is done based on the type of variants, since SNPRelate can handle different variants beyond bi-allelic SNPs. The filtered output is called an_ac_filtered.vcf. This script also creates the directories for the following steps and removes some of the files that are no longer needed to save space.
3. __Convert vcf -> gds for SNPRelate analysis__: The r script 03_convert_gds.R loads the SNPRelate package which reads the filtered .vcf file from step 2 and outputs it into SNPRelate's .gds format (into 03_gds_pruning dir).
4. __Run LD Pruning and PCA__: The r script 04_pruning_pca.R open the .gds file from step 3, performs linkage pruning to remove variants that are likely not inherented independendly (linkage disquilibrium) and runs the PCA on the remaning variants. Plots of the first 4 PCs and a Scree plot are printed into 04_pca/.
5. __LD decay__: To find an appropriate threshold for the LD pruning, a LD decay plot can be used. With the here used data and the SNPRelate package, there was no LD decay plot that looked as expected. See code to create the plot and further detail in 05_ld_decay.R.

### To Do & Known Issues (if there was more time)
- add Snakemake
- make pruning step a separated r script
- LD Pruning threshold as of now is not based on an LD decay plot of the used data
- redo the analysis with plink and try to solve the LD decay plot there => then run analysis with threshold that matches the plot
- sliding window PCA to look at Chromosome parts isolated and/or PCA per Chromosome instead of having one PCA for the whole genome

***************************************

# Data
The .vcf file contains variants called for 15 samples from 3 population, and one additonal sample ("Naxos2"), which is an outgroup for phylogenetics that we will not use for the PCA. The raw data is not provided here, reproduction is not possible without access to the university server.

```
less 01_raw/ProjTaxa.vcf
```

## Explaination of columns in "FORMAT" column:
```
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=SB,Number=4,Type=Integer,Description="Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias.">
```

GT is the genotype (reference or alternative), AD is the number of reads supporting each allele, DP total reads covering the site, GQ is the confidence of the genotype

## Explaination of columns in "INFO" column:
```
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
```

Present in our data are AC (Allele Count) and AN (Allele Number). AC gives the number of one genotype (for a SNP one Nucleotide) across the samples, that is not the reference. Since there can be multiple alterantives to the reference, the AC can give multiple values (matches order of alternatives). The Allele Number is the total number of how often that site was called (regardless of which allele it is, just how often we have data for the site). => The AN could be useful for filtering to exclude sites for which we do not have much data, and AC could be useful to filter very rare variants.

## Columns headers in the vcf file
```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  8N05240 8N05890 8N06612 8N73248 8N73604 K006    K010    K011    K015    K019    Lesina_280      Lesina_281      Lesina_282      Lesina_285      Lesina_286        Naxos2
```

Everything after FORMAT are the individual samples. The FILTER column was not used at all (just "."). The INFO columns gives the values AC and AN for each variant.

Further Info on VCF file structure: https://samtools.github.io/hts-specs/VCFv4.5.pdf

## Example row
```
# chr5    1002    .       T       C       1711.32 .       AC=2;AN=4       GT:AD:DP:GQ:PL  ./.:71,0:71:.:. ./.:120,0:120:.:.       ./.:94,0:94:.:. 0/0:72,0:72:0:0,0,25    ./.:124,0:124:.:.       ./.:119,0:119:.:.       ./.:86,0:86:.:.       ./.:72,0:72:.:. ./.:77,0:77:.:. ./.:109,0:109:.:.       ./.:29,0:29:.:. ./.:12,0:12:.:. 1/1:1,9:10:2:254,2,0    ./.:8,0:8:.:.   ./.:20,0:20:.:. ./.:23,0:23:.:.
```
