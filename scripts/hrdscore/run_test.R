#!/usr/bin/env Rscript

' run_test.R

Usage: run_test.R -l LOH -o OUTPUT -r GENOME [ --hrdtools HRDTOOLS ]

Options:
    -l --loh LOH            Path to input LOH segments file

    -o --output OUTPUT      Path to output (tsv file) summary of HRD scores

    -r --genome GENOME      Name of genome (i.e. hg19)

    --hrdtools HRDTOOLS     If provided, runs in debug mode, where HRDTOOLS is the path to
                            the HRDtools package. The package will be loaded in using 
                            devtools::load_all() instead of library(hrdtools).
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

print('loading packages...')

if (!is.null(args[['hrdtools']])) {
    library(devtools)
    load_all(args[['hrdtools']])
    n <- readLines(paste(args[['hrdtools']], 'NAMESPACE', sep = '/'))
    for (p in unique(gsub('.*?\\((.*?)[,)].*', '\\1', n[grepl('^import', n)]))) library(p, character.only = TRUE)
} else {
    suppressMessages(library(hrdtools))
}

print('done')

print('Running HRDtools test')

loh_ranges <- read_tsv(args[['loh']]) %>% as.data.frame %>% GRanges

run_test(
    loh_ranges, 
    loh.col = 'lohtype',
    cnv.col = 'copy_number',
    genome = args[['genome']]
) %>%
    as_tibble() %>%
    write_tsv(args[['output']])

print(paste0('Done. Wrote output to ', args[['output']]))
