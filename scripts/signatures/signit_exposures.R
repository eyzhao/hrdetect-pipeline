' signit_exposures.R

Computes mutation signatures using SignIT. Accepts a mutation catalog as input.
Mutation catalogs can be output from get_mutation_catalog.R.

Can be run with different reference mutation signature tables.
- wtsi_30_snv_signatures: based on COSMIC 30-signature reference
- wtsi_30_snv_signatures_exome: as above, but corrected for trinucleotide frequency of exomes
- nikzainal_sv_signatures: breast cancer SV signatures (Nik-Zainal et al. 2016, https://www.nature.com/articles/nature17676)

Usage: signit_exposures.R -c CATALOG -o OUTPUT [ -s SIGNIT -r REFERENCE ]

Options:
    -c --catalog CATALOG        Path to mutation catalog. TSV with two columns: mutation_type (str) and count (int).

    -o --output OUTPUT          Path to .Rds file containing the serialized exposures object

    -s --signit SIGNIT          Path to SignIT R library files. If not provided, will assume that SignIT is installed
                                and load the package using library(signit) instead.

    -r --reference REFERENCE    Name of the reference mutation signature dataset to use. The default is called
                                wtsi_30_snv_signatures. However, for SV signatures, one must use
                                nikzainal_sv_signatures instead.

Examples:
    To get mutation catalogs:
    Rscript signatures/get_mutation_catalog.R -v somatic_variants.vcf -c mutation_catalog.tsv

    To compute signatures:
    Rscript signatures/signit_exposures.R -c mutation_catalog.tsv -o signit_output.Rds

    Using a custom SignIT path:
    Rscript signatures/signit_exposures.R -c mutation_catalog.tsv -o signit_output.Rds -s path/to/signit/package

    Analysis of exomes:
    Rscript signatures/signit_exposures.R -c mutation_catalog.tsv -o signit_output.Rds -r wtsi_30_snv_signatures_exome

    Analysis of SV signatures:
    Rscript sv/calculate_sv_catalog.R -i sv_calls.tsv -o sv_catalog.tsv
    Rscript signatures/signit_exposures.R -c sv_catalog.tsv -o signit_output.Rds -r nikzainal_sv_signatures
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

exposures <- get_exposures(file = args[['catalog']], reference_signatures = reference, n_chains = 4, n_cores = 4)

print('SignIT Analysis Complete.')

saveRDS(exposures, file = args[['output']])

print(paste0('Results saved as RDS serialized object in ', args[['output']]))
