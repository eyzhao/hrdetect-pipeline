' hrdetect_input_table.R

Usage: hrdetect_input_table.R (-i PATH)... -o OUTPUT

Options:
    -i --input PATH         Include one of these per file to merge
    -o --output OUTPUT      Path to output file with merged data

Example:
    hrdetect_input_table.R -i path/to/hrdscore -i path/to/microhomology
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
