from glob import glob
import re
import pandas as pd
import itertools
import sys
import glob

### CONFIGS ###

PROJECT_DIR = '/projects/ezhao_prj/analyses/hrdetect-pipeline'
PURDOM_BOOTSTRAP_ITERATIONS = '1000'
SIGNIT_PATH = '/projects/ezhao_prj/papers/SignIT-paper/analysis/scripts/SignIT'
HRDTOOLS_PATH = '/projects/ezhao_prj/dependencies/packages/hrdtools'

newline_re = re.compile(r'\s*\n\s*')
def replace_newlines(input_string):
    return newline_re.sub(' ', input_string)

#########################################################
### Load Modules for Project-Specific File Structures ###
#########################################################

target_functions = {}

include: "projects/compass.smk"
include: "projects/pog.smk"

#######################
### Loading Targets ###
#######################

def glob_to_df(glob_string, regex, column_names):
    regex_obj = re.compile(regex)
    print('Running glob: ' + glob_string)
    glob_output = glob.iglob(glob_string)
    regex_searches = (regex_obj.search(s) for s in glob_output)
    regex_matches = [[r.group(i) for i in range(1, len(column_names)+1)] for r in regex_searches]
    df = pd.DataFrame(
        regex_matches
    )
    df = df.iloc[:, range(len(column_names))]
    df.columns = column_names
    print('Glob complete.')
    return df

def get_aggregation_targets(wildcards):
    target_loading_function = \
        target_functions[wildcards.project][wildcards.subproject][wildcards.signature_type]

    targets_df = target_loading_function(wildcards)

    targets = 'output/' + \
        wildcards.project + '/' + \
        wildcards.subproject + '/patients/' + \
        targets_df.patient + '/' + \
        targets_df.samples + '/' + \
        target_path_endings[wildcards.signature_type]

    return(targets)

target_path_endings = {
    'snv_signatures': 'snv_signatures/signit_hrd_signatures.tsv',
    'microhomology_scores': 'indel_microhomology/microhomology_scores.tsv',
    'sv_signatures': 'sv_signatures/sv_signit_hrd_signatures.tsv',
    'hrd_scores': 'hrd_scores/hrd_score.tsv'
}

##################################
### The Main HRDetect Pipeline ###
##################################

rule all:
    input:
        'output/compass/all/hrdetect/hrdetect_output_table.tsv',
        'output/pog/all/hrdetect/hrdetect_output_table.tsv',

rule aggregation_paths:
    input:
        get_aggregation_targets
    output:
        'output/{project}/{subproject}/hrdetect/{signature_type}_paths.tsv'
    log:
        'logs/{project}/{subproject}/hrdetect/{signature_type}_paths.log'
    run:
        f = open(output[0], 'w')
        f.write('\n'.join(input))
        f.close()

rule aggregate:
    input:
        'output/{project}/{subproject}/hrdetect/{signature_type}_paths.tsv'
    output:
        'output/{project}/{subproject}/hrdetect/{signature_type}_aggregated.tsv'
    log:
        'logs/{project}/{subproject}/hrdetect/{signature_type}_aggregated.log'
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/basic-functions/row_bind_tables.R -i {input} -o {output}'

###############################
### SNV Mutation Signatures ###
###############################

rule mutation_catalogs:
    input:
        'data/{project}/{subproject}/patients/{patient}/{sample}/somatic_snvs.vcf'
    output:
        'output/{project}/{subproject}/patients/{patient}/{sample}/snv_signatures/mutation_catalog.tsv',
    log:
        'logs/{project}/{subproject}/patients/{patient}/{sample}/snv_signatures/mutation_catalog.log'
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/signatures/get_mutation_catalog.R -v {input} -c {output} -r hg19'

rule snv_signatures_signit_exposures:
    input:
        'output/{project}/{subproject}/patients/{patient}/{sample}/snv_signatures/mutation_catalog.tsv',
    output:
        'output/{project}/{subproject}/patients/{patient}/{sample}/snv_signatures/signit_output.Rds',
    log:
        'logs/{project}/{subproject}/patients/{patient}/{sample}/snv_signatures/signit_output.log',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/signatures/signit_exposures.R -c {input} -o {output} -s ' + SIGNIT_PATH

rule snv_signatures_signit_summary:
    input:
        'output/{project}/{subproject}/patients/{patient}/{sample}/snv_signatures/signit_output.Rds',
    output:
        exposure='output/{project}/{subproject}/patients/{patient}/{sample}/snv_signatures/signit_exposure_summary.tsv',
        fraction='output/{project}/{subproject}/patients/{patient}/{sample}/snv_signatures/signit_exposure_summary_fractions.tsv',
    log:
        'logs/{project}/{subproject}/patients/{patient}/{sample}/snv_signatures/signit_exposure_summary.log',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/signatures/signit_summary_table.R -i {input} -o {output.exposure} --signit ' + SIGNIT_PATH + \
        ' && Rscript ' + PROJECT_DIR + '/scripts/signatures/signit_summary_table.R -i {input} -o {output.fraction} --fraction --signit ' + SIGNIT_PATH

rule snv_hrd_signatures:
    input:
        'output/{project}/{subproject}/patients/{patient}/{sample}/snv_signatures/signit_exposure_summary.tsv',
    output:
        'output/{project}/{subproject}/patients/{patient}/{sample}/snv_signatures/signit_hrd_signatures.tsv',
    log:
        'logs/{project}/{subproject}/patients/{patient}/{sample}/snv_signatures/signit_hrd_signatures.log',
    shell:
        replace_newlines(
            """ Rscript -e '
                    library(tidyverse);
                    read_tsv("{input}") %>%
                    filter(signature %in% c("Signature 3", "Signature 8")) %>%
                    select(signature, median_exposure) %>%
                    spread(signature, median_exposure) %>% write_tsv("{output}")
            '
            """.replace('\n', ' ')
        )

##############################
### SV Mutation Signatures ###
##############################

rule sv_catalogs:
    input:
        'data/{project}/{subproject}/patients/{patient}/{sample}/somatic_sv.tsv'
    output:
        'output/{project}/{subproject}/patients/{patient}/{sample}/sv_signatures/sv_catalog.tsv'
    log:
        'logs/{project}/{subproject}/patients/{patient}/{sample}/sv_signatures/sv_catalog.log'
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/sv/calculate_sv_catalog.R -i {input} -o {output}'

rule sv_signatures_signit_exposures:
    input:
        'output/{project}/{subproject}/patients/{patient}/{sample}/sv_signatures/sv_catalog.tsv',
    output:
        'output/{project}/{subproject}/patients/{patient}/{sample}/sv_signatures/sv_signit_output.Rds',
    log:
        'logs/{project}/{subproject}/patients/{patient}/{sample}/sv_signatures/sv_signit_output.log',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/signatures/signit_exposures.R -c {input} -o {output} -r nikzainal_sv_signatures -s ' + SIGNIT_PATH

rule sv_signatures_signit_summary:
    input:
        'output/{project}/{subproject}/patients/{patient}/{sample}/sv_signatures/sv_signit_output.Rds',
    output:
        exposure='output/{project}/{subproject}/patients/{patient}/{sample}/sv_signatures/sv_signit_exposure_summary.tsv',
        fraction='output/{project}/{subproject}/patients/{patient}/{sample}/sv_signatures/sv_signit_exposure_summary_fractions.tsv',
    log:
        'logs/{project}/{subproject}/patients/{patient}/{sample}/sv_signatures/sv_signit_exposure_summary.log',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/signatures/signit_summary_table.R -i {input} -o {output.exposure} --signit ' + SIGNIT_PATH + \
        ' && Rscript ' + PROJECT_DIR + '/scripts/signatures/signit_summary_table.R -i {input} -o {output.fraction} --fraction --signit ' + SIGNIT_PATH

rule sv_hrd_signatures:
    input:
        'output/{project}/{subproject}/patients/{patient}/{sample}/sv_signatures/sv_signit_exposure_summary.tsv',
    output:
        'output/{project}/{subproject}/patients/{patient}/{sample}/sv_signatures/sv_signit_hrd_signatures.tsv',
    log:
        'logs/{project}/{subproject}/patients/{patient}/{sample}/sv_signatures/sv_signit_hrd_signatures.log',
    shell:
        replace_newlines(
            """ Rscript -e '
                    library(tidyverse);
                    read_tsv("{input}") %>%
                    filter(signature %in% c("Rearrangement Signature 3", "Rearrangement Signature 5")) %>%
                    select(signature, median_exposure) %>%
                    spread(signature, median_exposure) %>% write_tsv("{output}")
            '
            """.replace('\n', ' ')
        )

###########################
### Indel Microhomology ###
###########################

rule indel_annotation:
    input:
        'data/{project}/{subproject}/patients/{patient}/{sample}/somatic_indels.vcf'
    output:
        'output/{project}/{subproject}/patients/{patient}/{sample}/indel_microhomology/indels_annotated.tsv',
    log:
        'logs/{project}/{subproject}/patients/{patient}/{sample}/indel_microhomology/indels_annotated.log',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/microhomology/flank_indels.R -v {input} -o {output}'

rule microhomology_calculation:
    input:
        'output/{project}/{subproject}/patients/{patient}/{sample}/indel_microhomology/indels_annotated.tsv',
    output:
        'output/{project}/{subproject}/patients/{patient}/{sample}/indel_microhomology/microhomology_scores.tsv',
    log:
        'logs/{project}/{subproject}/patients/{patient}/{sample}/indel_microhomology/microhomology_scores.log',
    shell:
        'Rscript scripts/microhomology/compute_microhomology_scores.R -i {input} -o {output}'

##################
### HRD Scores ###
##################

rule hrd_scores:
    input:
        'data/{project}/{subproject}/patients/{patient}/{sample}/segments.tsv'
    output:
        'output/{project}/{subproject}/patients/{patient}/{sample}/hrd_scores/hrd_score.tsv'
    log:
        'logs/{project}/{subproject}/patients/{patient}/{sample}/hrd_scores/hrd_score.log'
    shell:
        'Rscript scripts/hrdscore/run_test.R -l {input} -o {output} -r hg19 --hrdtools ' + HRDTOOLS_PATH

#######################
### HRDetect Scores ###
#######################

rule hrdetect_input_table:
    input:
        snv = 'output/{project}/{subproject}/hrdetect/snv_signatures_aggregated.tsv',
        sv = 'output/{project}/{subproject}/hrdetect/sv_signatures_aggregated.tsv',
        hrd = 'output/{project}/{subproject}/hrdetect/hrd_scores_aggregated.tsv',
        microhomology = 'output/{project}/{subproject}/hrdetect/microhomology_scores_aggregated.tsv',
    output:
        'output/{project}/{subproject}/hrdetect/hrdetect_input_table.tsv',
    log:
        'logs/{project}/{subproject}/hrdetect/hrdetect_input_table.log',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/hrdetect/hrdetect_input_table.R -i {input.snv} -i {input.sv} -i {input.hrd} -i {input.microhomology} -o {output}'

rule hrdetect_output_table:
    input:
        'output/{project}/{subproject}/hrdetect/hrdetect_input_table.tsv',
    output:
        'output/{project}/{subproject}/hrdetect/hrdetect_output_table.tsv',
    log:
        'output/{project}/{subproject}/hrdetect/hrdetect_output_table.tsv',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/hrdetect/compute_hrdetect.R -i {input} -o {output}'

#########################
### Report Generation ###
#########################
