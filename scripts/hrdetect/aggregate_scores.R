' aggregate_scores.R

Usage: aggregate_scores.R -p PATIENT -s SAMPLE --hrd HRD --snv SNV --sv SV --mh MH --out OUT

Options:
    -p --patient PATIENT    Unique patient ID
    -s --sample SAMPLE      Unique sample ID
    --hrd HRD               Path to HRD score output for a case
    --snv SNV               Path to SNV signatures output for a case
    --sv SV                 Path to SV signatures output for a case
    --mh MH                 Path to indel microhomology output for a case
    --out OUT               Path to output file, with data aggregated across these tables
' -> doc

library(docopt)
args <- docopt(doc)

library(tidyverse)

hrd_table <- read_tsv(args[['hrd']]) %>%
    select(`HRD score` = total) %>%
    gather(signature, value)

snv_table <- read_tsv(args[['snv']]) %>%
    select(signature, value = mean_exposure) %>%
    filter(signature %in% c('Signature 3', 'Signature 8'))

sv_table <- read_tsv(args[['sv']]) %>%
    select(signature, value = mean_exposure) %>%
    filter(signature %in% c('Rearrangement Signature 3', 'Rearrangement Signature 5'))

mh_table <- read_tsv(args[['mh']]) %>%
    select(`Deletion Microhomology Proportion` = deletion_microhomology_proportion) %>%
    gather(signature, value)

combined <- hrd_table %>%
    bind_rows(snv_table) %>%
    bind_rows(sv_table) %>%
    bind_rows(mh_table) %>%
    mutate(
        patient = args[['patient']],
        sample = args[['sample']]
    ) %>%
    select(patient, sample, signature, value)

combined %>%
    write_tsv(args[['out']])
