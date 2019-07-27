' signit_summary_table.R

Extracts data from serialized SignIT output file into a tabular format.
After running signit_exposures.R, you will obtain a .Rds file.
This script reads tabular output from that Rds file.

Usage: signit_summary_table.R -i INPUT -o OUTPUT [ -s SIGNIT --fraction ]

Options:
    -i --input INPUT            Path to serialized output from SignIT (.Rds file)

    -o --output OUTPUT          Path to output summary stats table (.tsv file)

    -s --signit SIGNIT          Path to SignIT R library files. If not provided, will assume that SignIT is installed
                                and load the package using library(signit) instead.

    --fraction                  If using this flag, summary table will report all values as exposure fractions instead
                                of number of mutations.

Examples:
    To get mutation catalogs:
    Rscript signatures/get_mutation_catalog.R -v somatic_variants.vcf -c mutation_catalog.tsv

    To compute signatures:
    Rscript signatures/signit_exposures.R -c mutation_catalog.tsv -o signit_output.Rds

    To extract tabular results:
    Rscript signatures/signit_summary_table.R -i signit_output.Rds -o signit_results.tsv
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
