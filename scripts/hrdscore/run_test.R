#!/usr/bin/env Rscript

' run_test.R

Runs HRD testing using the HRDtools package. Given a file with CNV/LOH segments, computes
HRD-LOH, HRD-TAI, and HRD-LST scores as described by Timms et al. (2014).
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4308910/.

The input LOH segments file must be structured calls with at least 5 columns.

- chr: The chromosome name
- start: Start position of CNV/LOH call
- end: End position of CNV/LOH call
- copy_number: The tumour copy number of the segment
- lohtype: The type of LOH state. Should be amongst the following:
    - ASCNA: Allele-specific copy number amplification
    - BCNA: Balanced copy number amplification
    - HET: Heterozygous (normal)
    - NLOH: Neutral LOH (loss of heterozygosity, but 2 copies present)
    - DLOH: Deletion LOH
    - ALOH: Amplification LOH

Usage: run_test.R -l LOH -o OUTPUT -r GENOME [ -h HRDTOOLS ]

Options:
    -l --loh LOH            Path to input LOH segments file

    -o --output OUTPUT      Path to output (tsv file) summary of HRD scores

    -r --genome GENOME      Name of genome (i.e. hg19)

    -h --hrdtools HRDTOOLS  If provided, runs in debug mode, where HRDTOOLS is the path to
                            the HRDtools package. The package will be loaded in using 
                            devtools::load_all() instead of library(hrdtools).

Examples:
    Rscript hrdscore/run_test.R -l loh_segs.tsv -o hrd_results.tsv -r hg19 --hrdtools path/to/hrdtools/package
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
