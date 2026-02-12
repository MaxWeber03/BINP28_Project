# Look at 01_raw/ProjTaxa.vcf
more 01_raw/ProjTaxa.vcf

# How many variants are present in the raw vcf?
# count lines that are not starting with #
grep -v ^\# 01_raw/ProjTaxa.vcf | wc -l
# 3816977

# How many chromosomes are present in the vcf?
# chromosomes start with "##contig=<ID=chr", so we can count those lines
grep ^\#\#contig\=\<ID\=chr* 01_raw/ProjTaxa.vcf | wc -l
# 30
