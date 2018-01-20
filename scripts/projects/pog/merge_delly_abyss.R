' merge_delly_abyss.R

Creates a merge between SV calls from Delly and SV calls from ABySS.

Usage: merge_delly_abyss.R -d DELLYPATH -a ABYSSPATH -o OUTPUT

Options:
    -d DELLYPATH        Path to Delly mutation calls.
                        These should be in VCF format and have a boolean
                        called Somatic which indicates whether the event
                        is somatic or not.

    -a ABYSSPATH        Path to ABySS mutation calls.

    -o OUTPUT           Path to output merged table.
' -> doc

library('docopt')
args <- docopt(doc)

library(tidyverse)
library(VariantAnnotation)

delly_path <- args[['d']]
abyss_path <- args[['a']]
output_path <- args[['o']]

print(paste('Merging SV calls from files:', delly_path, 'and', abyss_path))
vcf <- readVcf(delly_path, 'hg19')
table <- cbind(as.data.frame(rowRanges(vcf)), as.data.frame(info(vcf)))
delly_somatic <- table[table$Somatic, ]

abyss <- read_tsv(abyss_path)
abyss_tidy <- abyss %>% 
    separate(col = breakpoint, 
             into = c('chr1', 'pos1', 'chr2', 'pos2'), 
             sep = '[\\|\\:]') %>%
    dplyr::mutate(abyss_index = row_number())

match <- plyr::ddply(delly_somatic, c('seqnames', 'start'), function(z) {
    z <- z[1, ]
    m <- abyss_tidy[((abyss_tidy$chr1 == z$seqnames & abyss_tidy$chr2 == z$CHR2)
                    | (abyss_tidy$chr2 == z$seqnames & abyss_tidy$chr1 == z$CHR2))
                    & (abs(as.numeric(z$start) - as.numeric(abyss_tidy$pos1)) < 20 
                       | abs(as.numeric(z$start) - as.numeric(abyss_tidy$pos2))  < 20 )
                    & (abs(as.numeric(z$END) - as.numeric(abyss_tidy$pos2)) < 20 
                       | abs(as.numeric(z$END) - as.numeric(abyss_tidy$pos1)) < 20), ]
    if(dim(m)[1] > 0) {
        abyss_index = m[['abyss_index']][1]
    } else {
        abyss_index = NA
    }
    z <- z %>% mutate(abyss_index = abyss_index)
    return(z)
}) %>% 
    filter(!is.na(abyss_index)) %>%
    inner_join(abyss_tidy, by = 'abyss_index')

print(match %>% as_tibble)

match %>% dplyr::select(chr1, 
                 pos1, 
                 chr2, 
                 pos2,
                 type = SVTYPE
                 ) %>%
    write_tsv(path = output_path)
