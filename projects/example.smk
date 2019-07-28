def get_all_example_cases(wildcards):
    snv = get_example_snv_targets()
    sv = get_example_sv_targets()
    mh = get_example_microhomology_targets()
    hrd = get_example_hrd_score_targets()

    merged = snv.merge(
        sv, on=['patient', 'samples']
    ).merge(
        mh, on=['patient', 'samples']
    ).merge(
        hrd, on=['patient', 'samples']
    )

    return merged[['patient', 'samples']]

def get_example_snv_targets():
    return glob_to_df(
        'examples/example/patient*/sample*/snv.vcf',
        'examples/example/(.*?)/(.*?)/snv.vcf',
        ['patient', 'samples'],
    )

def get_example_microhomology_targets():
    return glob_to_df(
        'examples/example/patient*/sample*/indel.vcf',
        'examples/example/(.*?)/(.*?)/indel.vcf',
        ['patient', 'samples'],
    )

def get_example_hrd_score_targets():
    return glob_to_df(
        'examples/example/patient*/sample*/cnv.txt',
        'examples/example/(.*?)/(.*?)/cnv.txt',
        ['patient', 'samples'],
    )

def get_example_sv_targets():
    return glob_to_df(
        'examples/example/patient*/sample*/sv.tsv',
        'examples/example/(.*?)/(.*?)/sv.tsv',
        ['patient', 'samples'],
    )


target_functions['example'] = get_all_example_cases

rule example_snvs:
    input:
        ancient('examples/example/{patient}/{sample}/snv.vcf')
    output:
        'data/example/patients/{patient}/{sample}/somatic_snvs.vcf'
    shell:
        'cp -f {input} {output}'

rule example_indels:
    input:
        ancient('examples/example/{patient}/{sample}/indel.vcf')
    output:
        'data/example/patients/{patient}/{sample}/somatic_indels.vcf'
    shell:
        'cp -f {input} {output}'


rule example_segments:
    input:
        ancient('examples/example/{patient}/{sample}/cnv.txt')
    output:
        'data/example/patients/{patient}/{sample}/segments.tsv'
    run:
        data = pd.read_csv(
            input[0],
            sep='\t',
            names=['chr', 'start', 'end', 'width', 'nvar', 'copy_number', 'lohtype', 'major_allele', 'minor_allele']
        )
        data.to_csv(output[0], sep='\t', index=False)

rule example_sv:
    input:
        ancient('examples/example/{patient}/{sample}/sv.tsv')
    output:
        'data/example/patients/{patient}/{sample}/somatic_sv.tsv'
    shell:
        'cp -f {input} {output}'
