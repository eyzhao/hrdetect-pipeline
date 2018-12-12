' hrdetect_input_table.R

Gathers data necessary for HRDetect calculation into a single table.
The necessary inputs are 
    1. HRD scores: output by hrdscore/run_test.R
    2. SNV signatures: output by running signatures/signit_summary_table.R on SNV mutation catalog from signatures/get_mutation_catalog.R
    3. SV signatures: output by running signatures/signit_summary_table.R on SV mutation catalog from sv/calculate_sv_catalog.R
    4. Microhomology scores: output by running microhomology/compute_microhomology_scores.R

Usage: hrdetect_input_table.R (-i PATH)... -o OUTPUT

Options:
    -i --input PATH         Include one of these per file to merge
    -o --output OUTPUT      Path to output file with merged data

Example:
    hrdetect/hrdetect_input_table.R -i hrd_results.tsv -i microhomology_results.tsv -i snv_signatures.tsv -i sv_signatures.tsv -o hrdetect_input.tsv
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

input_paths <- tibble(path = args[['input']] %>% unique)

plyr::ddply(input_paths, 'path', function(p) {
    path = p$path
    read_tsv(path) %>%
        separate(paths, c('f1', 'f2', 'f3', 'f4', 'patient', 'sample', 'data_type', 'filename'), sep = '\\/') %>%
        select(-f1, -f2, -f3, -f4, -filename) %>%
        gather(variable, value, -patient, -sample, -data_type)
}) %>% 
    select(-path) %>%
    write_tsv(args[['output']])
