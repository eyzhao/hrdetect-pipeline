compass_project_root = '/projects/analysis/analysis30/BIOAPPS-7932'

def get_compass_snv_indel_targets(wildcards):
    return glob_to_df(
        '/projects/analysis/analysis30/BIOAPPS-7932/*/*.somatic.vcf',
        '/projects/analysis/analysis30/BIOAPPS-7932/(.*?)/(.*?).somatic.vcf',
        ['patient', 'samples']
    )

def get_compass_segment_targets(wildcards):
    return glob_to_df(
        '/projects/analysis/analysis30/BIOAPPS-7932/*/*.segments.txt',
        '/projects/analysis/analysis30/BIOAPPS-7932/(.*?)/(.*?).segments.txt',
        ['patient', 'samples']
    )

def get_compass_sv_targets(wildcards):
    return glob_to_df(
        '/projects/analysis/analysis30/BIOAPPS-7932/*/*.annotatedSV.tsv',
        '/projects/analysis/analysis30/BIOAPPS-7932/(.*?)/(.*?).annotatedSV.tsv',
        ['patient', 'samples']
    )

target_functions['compass'] = {
    'all': {
        'snv_signatures': get_compass_snv_indel_targets,
        'sv_signatures': get_compass_sv_targets,
        'hrd_scores': get_compass_segment_targets,
        'microhomology_scores': get_compass_snv_indel_targets,
    }
}

rule compass_snvs_and_indels:
    input:
        ancient(compass_project_root + '/{patient}/{sample}.somatic.vcf')
    output:
        snv='data/compass/all/patients/{patient}/{sample}/somatic_snvs.vcf',
        indel='data/compass/all/patients/{patient}/{sample}/somatic_indels.vcf',
    log:
        'logs/data/compass/all/patients/{patient}/{sample}/somatic_snvs.log',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/projects/compass/vcf_to_snv_indel.R -v {input} -s {output.snv} -i {output.indel}'

rule compass_segments:
    input:
        ancient(compass_project_root + '/{patient}/{sample}.segments.txt')
    output:
        'data/compass/all/patients/{patient}/{sample}/segments.tsv'
    log:
        'logs/data/compass/all/patients/{patient}/{sample}/segments.log'
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/projects/compass/format_segments.R -i {input} -o {output}'

rule compass_sv:
    input:
        ancient(compass_project_root + '/{patient}/{sample}.annotatedSV.tsv')
    output:
        'data/compass/all/patients/{patient}/{sample}/somatic_sv.tsv'
    log:
        'logs/data/compass/all/patients/{patient}/{sample}/somatic_sv.log'
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/projects/compass/format_sv.R -i {input} -o {output}'

