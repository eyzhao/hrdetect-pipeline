' vcf_to_snv_indel.R

Splits a VCF file into separate SNVs and Indel VCF files.

Usage: vcf_to_snv_indel.R -v VCF -s SNV -i INDEL

Options:
    -v --vcf VCF        Path to input VCF file
    -s --snv SNV        Path to output SNV VCF file
    -i --indel INDEL    Path to output Indel VCF file
' -> doc

library(docopt)
args <- docopt(doc)

library(VariantAnnotation)
library(tidyverse)

vcf <- readVcf(args[['vcf']])

snv <- vcf[isSNV(vcf, singleAltOnly=TRUE)]
indel <- vcf[!isSNV(vcf, singleAltOnly=FALSE)]

writeVcf(snv, args[['snv']])
print(paste0('Wrote SNVs to ', args[['snv']]))

writeVcf(indel, args[['indel']])
print(paste0('Wrote Indels to ', args[['indel']]))
