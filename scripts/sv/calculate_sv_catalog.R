' calculate_sv_catalog.R

Usage: calculate_sv_catalog.R -i INPUT -o OUTPUT

Options:
    -i --input INPUT        Table of SV data with at minimum the following columns
                                chr1: chromosome of the first breakpoint
                                chr2: chromosome of the second breakpoint
                                pos1: position of the first breakpoint
                                pos2: position of the second breakpoint
                                type: SV type with values DEL, DUP, INV, TRA
    
    -o --output OUTPUT      Path to output mutation catalog table
' -> doc

library(docopt)
args <- docopt(doc)

library('tidyverse')
library('copynumber')
library('GenomicRanges')

SV_BREAK_NAMES <- c('1-10kb', '10-100kb', '100kb-1Mb', '1Mb-10Mb', '>10Mb')
SV_LENGTH_BREAKS <- c(1000, 10000, 100000, 1000000, 10000000, Inf)

label_clustered <- function(sv_table) {
    id_table <- sv_table %>%
        mutate(id = 1:n())

    breakpoints <- bind_rows(
        id_table %>% 
            dplyr::select(id, chr = chr1, pos = pos1), 
        id_table %>% 
            dplyr::select(id, chr = chr2, pos = pos2)
    ) %>%
    mutate(
        chr = gsub('chr', '', chr), 
        chr = gsub('23', 'X', chr),
        chr = gsub('24', 'Y', chr),
        chr = factor(as.character(chr), levels=c(as.character(1:22), 'X', 'Y')),
        pos = as.numeric(pos)
    ) %>%
    arrange(chr, pos) %>%
    group_by(chr) %>% 
    mutate(
        prev_pos = c(0, pos[1:n()-1]), 
        distance = pos - prev_pos
    ) %>%
    ungroup()

    mean_distance <- mean(breakpoints$distance)

    print('clustered regions')

    clustered_regions <- breakpoints %>% 
        dplyr::select(chr, pos, distance) %>% 
        as.data.frame %>%
        pcf(kmin = 10, gamma = 25) %>%
        filter(mean < mean_distance / 10) %>%
        dplyr::select(chr = chrom, start = start.pos, end = end.pos) 

    if (dim(clustered_regions)[1] == 0) {
        breakpoints$clustered <- FALSE
    } else {
        clustered_breakpoint_indices <- breakpoints %>% 
            mutate(start = pos, end = pos + 1) %>% 
            dplyr::select(chr, start, end, id) %>% 
            GRanges() %>% 
            findOverlaps(clustered_regions %>% GRanges()) %>% 
            queryHits()

        breakpoints <- breakpoints %>% 
            mutate(
                clustered = if_else(
                    condition = row_number() %in% clustered_breakpoint_indices, 
                    true = TRUE,
                    false = FALSE
                )
            )
    }

    merged <- breakpoints %>% 
        group_by(id) %>% 
        summarise(clustered = any(clustered)) %>% 
        inner_join(id_table, by = 'id')

    return(merged)
}

calculate_sv_catalog <- function(sv_table) {
    base_table <- tibble(
        length = rep(c(rep(SV_BREAK_NAMES, 3), ''), 2) %>% as.character,
        type = rep(c(rep('DEL', 5), rep('DUP', 5), rep('INV', 5), 'TRA'), 2) %>% as.character,
        clustered = c(rep('clustered', 16), rep('non-clustered', 16)) %>% as.character
    )

    if (dim(sv_table)[1] <= 0) {
        sv_catalog <- base_table %>%
            mutate(
                count = 0
            )
    } else {
        sv_catalog <- sv_table %>% 
            label_clustered %>% 
            mutate(
                length = ifelse(type == 'TRA', 0, pos2-pos1), 
                length = cut(length, SV_LENGTH_BREAKS, labels = SV_BREAK_NAMES) %>% as.character
            ) %>% 
            replace_na(
                list(length = '')
            ) %>% 
            group_by(clustered, type, length) %>% 
            summarise(
                count = n()
            ) %>% 
            ungroup() %>% 
            mutate(
                clustered = ifelse(clustered, 'clustered', 'non-clustered')
            ) %>% 
            right_join(base_table, by = c('clustered', 'type', 'length')) %>% 
            mutate(
                count = ifelse(is.na(count), 0, count)
            )
    }

    print(sv_catalog %>% as.data.frame())
    return(sv_catalog %>%
        unite(
            'mutation_type',
            clustered, type, length,
            sep = '|'
        ))
}

sv_table <- read_tsv(args[['input']], col_types = cols(
        chr1 = col_character(),
        pos1 = col_integer(),
        chr2 = col_character(),
        pos2 = col_integer(),
        type = col_character()
    )) %>%
    calculate_sv_catalog() %>%
    write_tsv(args[['output']])

message(paste0('Wrote SV catalog to ', args[['output']]))
