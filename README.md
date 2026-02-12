Max Weber BINP28 M. Sc. Bioinformatics Programme Lund University, Feb 2026


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

## Columns headers in the vcf file
```#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  8N05240 8N05890 8N06612 8N73248 8N73604 K006    K010    K011    K015    K019    Lesina_280      Lesina_281      Lesina_282      Lesina_285      Lesina_286        Naxos2```

Everything after FORMAT are the individual samples.

## Example row
```# chr5    1002    .       T       C       1711.32 .       AC=2;AN=4       GT:AD:DP:GQ:PL  ./.:71,0:71:.:. ./.:120,0:120:.:.       ./.:94,0:94:.:. 0/0:72,0:72:0:0,0,25    ./.:124,0:124:.:.       ./.:119,0:119:.:.       ./.:86,0:86:.:.       ./.:72,0:72:.:. ./.:77,0:77:.:. ./.:109,0:109:.:.       ./.:29,0:29:.:. ./.:12,0:12:.:. 1/1:1,9:10:2:254,2,0    ./.:8,0:8:.:.   ./.:20,0:20:.:. ./.:23,0:23:.:.```

**********************************************

# Notes
- Try to look at other papers to find tools
- the SNPs will probably need to be heavily filtered
- PCA very much depends on input data quality

## Worflow outline
1. Filter VCF for only good SNPs => looks at papers for threshold and tools
2. Run linkage pruning => SNPs should be inherented intendepently. Looking at linkage disequilibrium would be good, but available time does not allow that.
3. Run PCA => find SNPs that are explaining variance between the samples/clusters
    - Expectation: samples should cluster together withhin each population

## To Do
- find software to filter SNPs, find approiate thresholds, remove Naxos2 outgroup


