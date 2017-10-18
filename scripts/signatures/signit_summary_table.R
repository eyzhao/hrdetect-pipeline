' signit_summary_table.R

Usage: signit_summary_table.R -i INPUT -o OUTPUT [ -s SIGNIT --fraction ]

Options:
    -i --input INPUT            Path to serialized output from SignIT (.Rds file)

    -o --output OUTPUT          Path to output summary stats table (.tsv file)

    -s --signit SIGNIT          Path to SignIT R library files. If not provided, will assume that SignIT is installed
                                and load the package using library(signit) instead.

    --fraction                  If using this flag, summary table will report all values as exposure fractions instead
                                of number of mutations.
' -> doc

library(docopt)
args <- docopt(doc)

if (is.null(args[['signit']])) {
    library(signit)
} else {
    library(devtools)
    library(tidyverse)
    library(rjags)
    library(nnls)
    library(dbscan)
    library(Rtsne)

    load_all(args[['signit']])
}

message('Reading SignIT data')

exposures <- readRDS(args[['input']])
summary_table <- get_exposure_summary_table(exposures, alpha = c(0, 0.05, 0.5), fraction = args[['fraction']])

message('Created summary table.')

summary_table %>% write_tsv(args[['output']])

print(paste0('Results saved as TSV at ', args[['output']]))
