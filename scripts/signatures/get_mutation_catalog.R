' run_deconstructsigs.R - Wrapper for deconstructsigs to compute mutation signatures

Usage: deconstructsigs_script.R ( -m MAF | -v VCF ) -c CATALOG [ -x -r REFERENCE ]

Options:
    -m MAF              SNV input in MAF format
    -v VCF              SNV input in VCF format
    -c CATALOG          Output to mutation catalog, computed by deconstructSigs
    -r REFERENCE        Reference genome name. Either hg19 or GRCh38. Default is hg19.

    -x --exome          Flag to indicate the variants were called from exome, which will
                            trigger the script to apply trinucleotide frequency correction
                            and adjust for mutation
' -> doc

library(docopt)
args <- docopt(doc)
print(args)

library(tidyverse)
library(deconstructSigs)
library(stringr)

if (! is.null(args[['m']])) {
    mutations <- read_tsv(args[['m']], 
                          col_types = cols(Chromosome = col_character(),
                                           Start_Position = col_number(),
                                           Reference_Allele = col_character(),
                                           Allele = col_character())) %>%
        select(chr = Chromosome,
               pos = Start_Position,
               ref = Reference_Allele,
               alt = Allele,
               sample = Tumor_Sample_Barcode) %>%
        mutate(chr = gsub('chr', '', chr)) %>%
        as.data.frame
} else {
    library(VariantAnnotation)
    vcf <- readVcf(args[['v']])
    mutations <- rowRanges(vcf)
    mutations <- mutations[
        elementNROWS(mutations$REF) == 1 &
        elementNROWS(mutations$ALT) == 1
    ]
    mutations$REF <- as.character(mutations$REF)
    mutations$ALT <- unlist(CharacterList(mutations$ALT))
    mutations <- as.data.frame(mutations) %>%
        dplyr::mutate(sample = args[['v']]) %>%
        dplyr::select(chr = seqnames,
               pos = start,
               ref = REF,
               alt = ALT,
               sample = sample) %>%
        as.data.frame()
}

print(mutations %>% as_tibble)

if (!is.null(args[['r']])) {
    if (args[['r']] == 'GRCh38') {
        library(BSgenome.Hsapiens.NCBI.GRCh38)

        catalog <- mut.to.sigs.input(mut.ref = mutations,
                                    sample.id = 'sample',
                                    chr = 'chr',
                                    pos = 'pos',
                                    ref = 'ref',
                                    alt = 'alt',
                                    bsg = BSgenome.Hsapiens.NCBI.GRCh38)
    } else {
        catalog <- mut.to.sigs.input(mut.ref = mutations,
                                    sample.id = 'sample',
                                    chr = 'chr',
                                    pos = 'pos',
                                    ref = 'ref',
                                    alt = 'alt')
    }
}

if (args[['exome']]) {
    correction = 'exome2genome'
    print('Running mutation signature detection using Exome to Genome trinucleotide frequency correction')
} else {
    correction = 'default'
}

catalog %>% t %>% as.data.frame %>% `colnames<-`(c('count')) %>% rownames_to_column('mutation_type') %>% write_tsv(args[['c']])

print(paste('Wrote catalogs to', args[['c']]))
