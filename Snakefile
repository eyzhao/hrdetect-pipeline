from glob import glob
import re
import pandas as pd
import itertools
import sys
import glob

### CONFIGS ###

PURDOM_BOOTSTRAP_ITERATIONS = '1000'
SIGNIT_PATH = 'git/SignIT'
HRDTOOLS_PATH = 'git/hrdtools'

newline_re = re.compile(r'\s*\n\s*')
def replace_newlines(input_string):
    return newline_re.sub(' ', input_string)

#########################################################
### Load Modules for Project-Specific File Structures ###
#########################################################

target_functions = {}

include: "projects/example.smk"
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

def get_hrdetect_targets(wildcards):
    print(wildcards)
    target_loading_function = \
        target_functions[wildcards.project]

    targets_df = target_loading_function(wildcards)
    print(targets_df)

    print('blah' + targets_df.patient + targets_df.samples)
    targets = 'output/' + \
        wildcards.project + '/' + \
        'patients/' + \
        targets_df.patient + '/' + \
        targets_df.samples + '/' + \
        'hrdetect_components.tsv'
    print(targets)

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
        'output/example/cohorts/all/hrdetect/hrdetect_output_table.tsv'

###############################
### SNV Mutation Signatures ###
###############################

rule mutation_catalogs:
    input:
        'data/{project}/patients/{patient}/{sample}/somatic_snvs.vcf'
    output:
        'output/{project}/patients/{patient}/{sample}/snv_signatures/mutation_catalog.tsv',
    shell:
        'Rscript scripts/signatures/get_mutation_catalog.R -v {input} -c {output} --ref hg19'

rule snv_signatures_signit_exposures:
    input:
        'output/{project}/patients/{patient}/{sample}/snv_signatures/mutation_catalog.tsv',
    output:
        'output/{project}/patients/{patient}/{sample}/snv_signatures/signit_output.Rds',
    shell:
        'Rscript scripts/signatures/signit_exposures.R -c {input} -o {output} -s ' + SIGNIT_PATH

rule snv_signatures_signit_summary:
    input:
        'output/{project}/patients/{patient}/{sample}/snv_signatures/signit_output.Rds',
    output:
        exposure='output/{project}/patients/{patient}/{sample}/snv_signatures/signit_exposure_summary.tsv',
        fraction='output/{project}/patients/{patient}/{sample}/snv_signatures/signit_exposure_summary_fractions.tsv',
    shell:
        'Rscript scripts/signatures/signit_summary_table.R -i {input} -o {output.exposure} --signit ' + SIGNIT_PATH + \
        ' && Rscript scripts/signatures/signit_summary_table.R -i {input} -o {output.fraction} --fraction --signit ' + SIGNIT_PATH

rule snv_hrd_signatures:
    input:
        'output/{project}/patients/{patient}/{sample}/snv_signatures/signit_exposure_summary.tsv',
    output:
        'output/{project}/patients/{patient}/{sample}/snv_signatures/signit_hrd_signatures.tsv',
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
        'data/{project}/patients/{patient}/{sample}/somatic_sv.tsv'
    output:
        'output/{project}/patients/{patient}/{sample}/sv_signatures/sv_catalog.tsv'
    shell:
        'Rscript scripts/sv/calculate_sv_catalog.R -i {input} -o {output}'

rule sv_signatures_signit_exposures:
    input:
        'output/{project}/patients/{patient}/{sample}/sv_signatures/sv_catalog.tsv',
    output:
        'output/{project}/patients/{patient}/{sample}/sv_signatures/sv_signit_output.Rds',
    shell:
        'Rscript scripts/signatures/signit_exposures.R -c {input} -o {output} --reference nikzainal_sv_signatures -s ' + SIGNIT_PATH

rule sv_signatures_signit_summary:
    input:
        'output/{project}/patients/{patient}/{sample}/sv_signatures/sv_signit_output.Rds',
    output:
        exposure='output/{project}/patients/{patient}/{sample}/sv_signatures/sv_signit_exposure_summary.tsv',
        fraction='output/{project}/patients/{patient}/{sample}/sv_signatures/sv_signit_exposure_summary_fractions.tsv',
    shell:
        'Rscript scripts/signatures/signit_summary_table.R -i {input} -o {output.exposure} --signit ' + SIGNIT_PATH + \
        ' && Rscript scripts/signatures/signit_summary_table.R -i {input} -o {output.fraction} --fraction --signit ' + SIGNIT_PATH

rule sv_hrd_signatures:
    input:
        'output/{project}/patients/{patient}/{sample}/sv_signatures/sv_signit_exposure_summary.tsv',
    output:
        'output/{project}/patients/{patient}/{sample}/sv_signatures/sv_signit_hrd_signatures.tsv',
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
        'data/{project}/patients/{patient}/{sample}/somatic_indels.vcf'
    output:
        'output/{project}/patients/{patient}/{sample}/indel_microhomology/indels_annotated.tsv',
    shell:
        'Rscript scripts/microhomology/flank_indels.R -v {input} -o {output}'

rule microhomology_calculation:
    input:
        'output/{project}/patients/{patient}/{sample}/indel_microhomology/indels_annotated.tsv',
    output:
        'output/{project}/patients/{patient}/{sample}/indel_microhomology/microhomology_scores.tsv',
    shell:
        'Rscript scripts/microhomology/compute_microhomology_scores.R -i {input} -o {output}'

##################
### HRD Scores ###
##################

rule hrd_scores:
    input:
        'data/{project}/patients/{patient}/{sample}/segments.tsv'
    output:
        'output/{project}/patients/{patient}/{sample}/hrd_scores/hrd_score.tsv'
    shell:
        'Rscript scripts/hrdscore/run_test.R -l {input} -o {output} --genome hg19 --hrdtools=' + HRDTOOLS_PATH

#######################
### HRDetect Scores ###
#######################

rule combine_scores:
    input:
        sv = 'output/{project}/patients/{patient}/{sample}/sv_signatures/sv_signit_exposure_summary.tsv',
        mh = 'output/{project}/patients/{patient}/{sample}/indel_microhomology/microhomology_scores.tsv',
        snv = 'output/{project}/patients/{patient}/{sample}/snv_signatures/signit_exposure_summary.tsv',
        hrd = 'output/{project}/patients/{patient}/{sample}/hrd_scores/hrd_score.tsv'
    output:
        'output/{project}/patients/{patient}/{sample}/hrdetect_components.tsv'
    shell:
        replace_newlines('''
            Rscript scripts/hrdetect/aggregate_scores.R
                -p {wildcards.patient}
                -s {wildcards.sample}
                --hrd {input.hrd}
                --snv {input.snv}
                --sv {input.sv}
                --mh {input.mh}
                --out {output}
        ''')

rule hrdetect_component_paths:
    input:
        get_hrdetect_targets
    output:
        'output/{project}/cohorts/{cohort}/hrdetect/hrdetect_input_table_paths.tsv',
    run:
        out = open(output[0], 'w')
        out.write('\n'.join(input))
        out.close()

rule aggregate_hrdetect_components:
    input:
        'output/{project}/cohorts/{cohort}/hrdetect/hrdetect_input_table_paths.tsv',
    output:
        'output/{project}/cohorts/{cohort}/hrdetect/hrdetect_input_table.tsv'
    shell:
        'Rscript scripts/basic-functions/row_bind_tables.R -i {input} -o {output}'

rule hrdetect_output_table:
    input:
        'output/{project}/cohorts/{cohort}/hrdetect/hrdetect_input_table.tsv',
    output:
        'output/{project}/cohorts/{cohort}/hrdetect/hrdetect_output_table.tsv',
    log:
        'output/{project}/cohorts/{cohort}/hrdetect/hrdetect_output_table.tsv',
    shell:
        'Rscript scripts/hrdetect/compute_hrdetect.R -i {input} -o {output}'

