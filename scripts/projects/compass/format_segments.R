' format_segments.R

Formats segment files from COMPASS project into the a TSV containing
columns chr, start, end, lohtype, copy_number, ready for import into
HRDtools.

Usage: format_segments.R -i INPUT -o OUTPUT

Options:
    -i --input INPUT        Path to input segment file from COMPASS.
    -o --output OUTPUT      Path to output segment file for HRDtools.
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

read_tsv(
    args[['input']],
    col_types = cols(labels = col_character())
) %>%
    select(
        chr = chrom,
        start = start.pos,
        end = end.pos,
        labels
    ) %>%
    separate(
        labels,
        c('minor', 'major'),
        sep = '\\.'
    ) %>%
    mutate(
        chr = gsub('chr', '', chr),
        chr = gsub('23', 'X', chr),
        chr = gsub('24', 'Y', chr),
        minor = as.numeric(minor),
        major = as.numeric(major)
    ) %>%
    filter(
        !is.na(minor) & !is.na(major)
    ) %>% 
    mutate(
        copy_number = minor + major,
        lohtype = case_when(
            minor == 0 & major == 0             ~ 'HOMD',
            minor == 0 & major == 1             ~ 'DLOH',
            minor == 0 & major == 2             ~ 'NLOH',
            minor == 0 & major > 2              ~ 'ALOH',
            minor == 1 & major == 1             ~ 'HET',
            minor != major                      ~ 'ASCNA',
            minor == major                      ~ 'BCNA'
        )
    ) %>%
    select(
        chr, start, end, copy_number, lohtype
    ) %>%
    write_tsv(args[['output']])
