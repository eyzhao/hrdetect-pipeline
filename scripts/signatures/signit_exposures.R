' signit_exposures.R

Usage: signit_exposures.R -c CATALOG -o OUTPUT [ -s SIGNIT -r REFERENCE ]

Options:
    -c --catalog CATALOG        Path to mutation catalog. TSV with two columns: mutation_type (str) and count (int).

    -o --output OUTPUT          Path to .Rds file containing the serialized exposures object

    -s --signit SIGNIT          Path to SignIT R library files. If not provided, will assume that SignIT is installed
                                and load the package using library(signit) instead.

    -r --reference REFERENCE    Name of the reference mutation signature dataset to use. The default is called
                                wtsi_30_snv_signatures. However, for SV signatures, one must use
                                nikzainal_sv_signatures instead.
' -> doc

library(docopt)
args <- docopt(doc)

if (is.null(args[['signit']])) {
    library(signit)
} else {
    library(devtools)
    library(tidyverse)
    library(rstan)
    library(nnls)
    library(dbscan)
    library(Rtsne)

    load_all(args[['signit']])
}

if (is.null(args[['reference']])) {
    ref_name = 'wtsi_30_snv_signatures'
} else {
    ref_name = args[['reference']]
}

reference <- get(data(list = ref_name))
print(reference)

catalog <- read_tsv(args[['catalog']])
exposures <- get_exposures(catalog, reference_signatures = reference, n_chains = 4, n_cores = 4)

print('SignIT Analysis Complete.')

saveRDS(exposures, file = args[['output']])

print(paste0('Results saved as RDS serialized object in ', args[['output']]))
