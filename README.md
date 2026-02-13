Max Weber BINP28 M. Sc. Bioinformatics Programme Lund University, Feb 2026

# Git Repository
Git Repository for this project can be found here: https://github.com/MaxWeber03/BINP28_Project.

I know that having both the report (latex code) and the project code in the same repo is not a clean solution, but I want to avoid bloating my github with multiple repos for a small project.

# Instructions Summary
We have a .vcf file without much futher detail on the data or the species. I am tasked with "Genetic relationship/Multivariate relationship from a PCA". We do not have much information. 3.5 days are available including writing of a two pages report and a two slide/3 min presentation.

# Data
The .vcf file contains variants called for 15 samples from 3 population, and one additonal sample ("Naxos2"), which is an outgroup for phylogenetics that we will not use for the PCA.

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

**********************************************

# Notes
- Try to look at other papers to find tools
- the SNPs will probably need to be heavily filtered
- PCA very much depends on input data quality

## Worflow outline
1. Filter VCF for only good SNPs => looks at papers and documentation for threshold
    - use bcftools, not vcftools since that is deprecated
    - https://samtools.github.io/bcftools/howtos/index.html
2. Run linkage pruning => SNPs should be inherented intendepently. Looking at linkage disequilibrium would be good, but available time does not allow that. (Will be done in SNPRelate in R)
3. Run PCA => find SNPs that are explaining variance between the samples/clusters
    - Expectation: samples should cluster together withhin each population
    - SNPRelate: https://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html


### Filtering/Setup
- Take away outgroup first!
- Quality Scores => look for threshold recommendations in bcftools
- Number of Reads => too many are artifacts, too few are untrsuted
    - look at distribution/histrogram of my data
        - Software for histrogram? => look for tools?
    - apply it as overall (per SNP)
        - not more then 2x mean/median is a usable value
- Only Bi-allelic SNPs or bi-allelic SNPs, multi-allelic SNPs, indels and structural variants?
    - focus on getting one done, we can compare both => go first with Bi-allelic
    - SNPRelate can handle INDELS and structural variants
        - I could try both if I have the time
        - no INDEL relalignment, that was done before variant calling
- How many PCs?
    - can be done through an analysis, but uncommen in this use case, number is going to be very high since we have so many dimensions => expect only 1-3% of variance explained per PC => hence the number of PCs does not have a huge effect on the results
- Minor Allel Frequency => don't worry about it for now
    - could be e.g. to give SNP confidence for downstream analysis

### VCF Filtering Details

- AC/AN => needs to be recalculated after removing the outgroup
    - bcftools filltags => find in docs
- Simple Options: Quality of SNP and of GT

### Further Ideas
- add Snakemake (if the available time allows)
