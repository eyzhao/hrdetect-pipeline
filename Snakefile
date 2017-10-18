from glob import glob
import re
import pandas as pd

### CONFIGS ###

PROJECT_DIR = '/projects/ezhao_prj/analyses/hrdetect-pipeline'
PURDOM_BOOTSTRAP_ITERATIONS = '1000'
SIGNIT_PATH = '/projects/ezhao_prj/papers/SignIT-paper/analysis/scripts/SignIT'
HRDTOOLS_PATH = '/projects/ezhao_prj/dependencies/packages/hrdtools'

#######################
### Loading Targets ###
#######################


#########################################################
### Load Modules for Project-Specific File Structures ###
#########################################################

include: "projects/compass.smk"


##################################
### The Main HRDetect Pipeline ###
##################################

#rule all:
#    input:
#
###############################
### SNV Mutation Signatures ###
###############################

rule mutation_catalogs:
    input:
        'data/{project}/{unique_sample_id}/somatic_snvs.vcf'
    output:
        'output/{project}/{unique_sample_id}/snv_signatures/mutation_catalog.tsv',
    log:
        'logs/{project}/{unique_sample_id}/snv_signatures/mutation_catalog.log'
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/signatures/get_mutation_catalog.R -v {input} -c {output} -r hg19'

rule snv_signatures_signit_exposures:
    input:
        'output/{project}/{unique_sample_id}/snv_signatures/mutation_catalog.tsv',
    output:
        'output/{project}/{unique_sample_id}/snv_signatures/signit_output.Rds',
    log:
        'logs/{project}/{unique_sample_id}/snv_signatures/signit_output.log',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/signatures/signit_exposures.R -c {input} -o {output} -s ' + SIGNIT_PATH

rule snv_signatures_signit_summary:
    input:
        'output/{project}/{unique_sample_id}/snv_signatures/signit_output.Rds',
    output:
        exposure='output/{project}/{unique_sample_id}/snv_signatures/signit_exposure_summary.tsv',
        fraction='output/{project}/{unique_sample_id}/snv_signatures/signit_exposure_summary_fractions.tsv',
    log:
        'logs/{project}/{unique_sample_id}/snv_signatures/signit_exposure_summary.log',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/signatures/signit_summary_table.R -i {input} -o {output.exposure} --signit ' + SIGNIT_PATH + \
        ' && Rscript ' + PROJECT_DIR + '/scripts/signatures/signit_summary_table.R -i {input} -o {output.fraction} --fraction --signit ' + SIGNIT_PATH


##############################
### SV Mutation Signatures ###
##############################

rule sv_catalogs:
    input:
        'data/{project}/{unique_sample_id}/somatic_sv.tsv'
    output:
        'output/{project}/{unique_sample_id}/sv_signatures/sv_catalog.tsv'
    log:
        'logs/{project}/{unique_sample_id}/sv_signatures/sv_catalog.log'
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/sv/calculate_sv_catalog.R -i {input} -o {output}'

rule sv_signatures_signit_exposures:
    input:
        'output/{project}/{unique_sample_id}/sv_signatures/sv_catalog.tsv',
    output:
        'output/{project}/{unique_sample_id}/sv_signatures/sv_signit_output.Rds',
    log:
        'logs/{project}/{unique_sample_id}/sv_signatures/sv_signit_output.log',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/signatures/signit_exposures.R -c {input} -o {output} -r nikzainal_sv_signatures -s ' + SIGNIT_PATH

rule sv_signatures_signit_summary:
    input:
        'output/{project}/{unique_sample_id}/sv_signatures/sv_signit_output.Rds',
    output:
        exposure='output/{project}/{unique_sample_id}/sv_signatures/sv_signit_exposure_summary.tsv',
        fraction='output/{project}/{unique_sample_id}/sv_signatures/sv_signit_exposure_summary_fractions.tsv',
    log:
        'logs/{project}/{unique_sample_id}/sv_signatures/sv_signit_exposure_summary.log',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/signatures/signit_summary_table.R -i {input} -o {output.exposure} --signit ' + SIGNIT_PATH + \
        ' && Rscript ' + PROJECT_DIR + '/scripts/signatures/signit_summary_table.R -i {input} -o {output.fraction} --fraction --signit ' + SIGNIT_PATH


###############################
### POG Indel Microhomology ###
###############################

rule indel_annotation:
    input:
        'data/{project}/{unique_sample_id}/somatic_indels.vcf'
    output:
        'output/{project}/{unique_sample_id}/indel_microhomology/indels_annotated.tsv',
    log:
        'logs/{project}/{unique_sample_id}/indel_microhomology/indels_annotated.log',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/microhomology/flank_indels.R -v {input} -o {output}'

rule microhomology_calculation:
    input:
        'output/{project}/{unique_sample_id}/indel_microhomology/indels_annotated.tsv',
    output:
        'output/{project}/{unique_sample_id}/indel_microhomology/microhomology_scores.tsv',
    log:
        'logs/{project}/{unique_sample_id}/indel_microhomology/microhomology_scores.log',
    shell:
        'Rscript scripts/microhomology/compute_microhomology_scores.R -i {input} -o {output}'

######################
### POG HRD Scores ###
######################

rule pog_hrd_scores:
    input:
        'data/{project}/{unique_sample_id}/segments.tsv'
    output:
        'output/{project}/{unique_sample_id}/hrd_scores/hrd_score.tsv'
    log:
        'logs/{project}/{unique_sample_id}/hrd_scores/hrd_score.log'
    shell:
        'Rscript scripts/hrdscore/run_test.R -l {input} -o {output} -r hg19 --hrdtools ' + HRDTOOLS_PATH

###########################
### POG HRDetect Scores ###
###########################



#########################
### Report Generation ###
#########################
