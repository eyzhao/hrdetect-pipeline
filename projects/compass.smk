compass_project_root = '/projects/analysis/analysis30/BIOAPPS-7932'

rule compass_snvs_and_indels:
    input:
        ancient(compass_project_root + '/{patient}/{sample}.somatic.vcf')
    output:
        snv='data/compass/{patient}_sample_{sample}/somatic_snvs.vcf',
        indel='data/compass/{patient}_sample_{sample}/somatic_indels.vcf',
    log:
        'logs/data/compass/{patient}_sample_{sample}/somatic_snvs.log',
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/projects/compass/vcf_to_snv_indel.R -v {input} -s {output.snv} -i {output.indel}'

rule compass_segments:
    input:
        ancient(compass_project_root + '/{patient}/{sample}.segments.txt')
    output:
        'data/{project}/{patient}_sample_{sample}/segments.tsv'
    log:
        'logs/data/{project}/{patient}_sample_{sample}/segments.log'
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/projects/compass/format_segments.R -i {input} -o {output}'

rule compass_sv:
    input:
        ancient(compass_project_root + '/{patient}/{sample}.annotatedSV.tsv')
    output:
        'data/{project}/{patient}_sample_{sample}/somatic_sv.tsv'
    log:
        'logs/data/{project}/{patient}_sample_{sample}/somatic_sv.log'
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/projects/compass/format_sv.R -i {input} -o {output}'



