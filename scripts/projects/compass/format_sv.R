' format_sv.R

Usage: format_sv.R -i INPUT -o OUTPUT

Options:
    -i --input INPUT        Path to COMPASS SV file
    -o --output OUTPUT      Path to output file for mutation signature analysis
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

read_tsv(args[['input']]) %>%
    filter(
        ! is.na(delly_id),
        ! is.na(crest_soft_clip1)
    ) %>%
    select(
        chr1 = chrom1, 
        pos1, 
        chr2 = chrom2, 
        pos2, 
        type
    ) %>%
    write_tsv(args[['output']])
