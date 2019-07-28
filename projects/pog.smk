pog_project_root = '/projects/POG/POG_data'

pog_custom_cohorts = pd.read_csv('config/pog_custom_cohorts.tsv', sep='\t')

class FileExistenceError(Exception):
    pass

class CohortError(Exception):
    pass

def is_integer(string):
    try:
        int(string)
        return True
    except ValueError:
        return False

def get_all_pogs():
    snv = get_pog_snv_targets()
    sv = get_pog_sv_targets()
    mh = get_pog_microhomology_targets()
    hrd = get_pog_hrd_score_targets()

    merged = snv.merge(
        sv, on=['patient', 'samples']
    ).merge(
        mh, on=['patient', 'samples']
    ).merge(
        hrd, on=['patient', 'samples']
    )

    return merged[['patient', 'samples']]

def get_pog_custom_cohort(cohort):
    return(
        pog_custom_cohorts[['patient', 'samples']][
            pog_custom_cohorts['cohort'] == cohort
        ]
    )

def get_pog_targets(wildcards):
    if wildcards.cohort == 'all':
        return(get_all_pogs())
    elif wildcards.cohort in pog_custom_cohorts.cohort.tolist():
        return(get_pog_custom_cohort(wildcards.cohort))
    elif wildcards.cohort in cancertype_cohorts:
        return(get_pog_cancertype_cohort(wildcards.cohort))
    else:
        raise(CohortError('Undefined cohort: {}'.format(wildcards.cohort)))

def get_pog_snv_targets():
    return glob_to_df(
        '/projects/POG/POG_data/POG*/wgs/*/strelka/bwa/results/passed.somatic.snvs.vcf',
        '/projects/POG/POG_data/(.*?)/wgs/(.*?)/strelka/bwa/results/passed.somatic.snvs.vcf',
        ['patient', 'samples'],
    )

def get_pog_microhomology_targets():
    return glob_to_df(
        '/projects/POG/POG_data/POG*/wgs/*/strelka/bwa/results/passed.somatic.indels.vcf',
        '/projects/POG/POG_data/(.*?)/wgs/(.*?)/strelka/bwa/results/passed.somatic.indels.vcf',
        ['patient', 'samples'],
    )

def get_pog_hrd_score_targets():
    return glob_to_df(
        '/projects/POG/POG_data/POG*/wgs/*/reviewed/loh/results/apolloh_out_segs.txt',
        '/projects/POG/POG_data/(.*?)/wgs/(.*?)/reviewed/loh/results/apolloh_out_segs.txt',
        ['patient', 'samples'],
    )

def get_pog_sv_targets():
    full_df = glob_to_df(
        '/projects/POG/POG_data/POG*/wgs/*_t_*_*_n_*',
        '/projects/POG/POG_data/(.*?)/wgs/(.*)',
        ['patient', 'samples'],
    )
    full_df[['t_prefix', 't', 't_lib', 'n_prefix', 'n', 'n_lib']] = full_df.samples.str.split('_', expand=True)

    delly_df = glob_to_df(
        '/projects/POG/POG_data/POG*/sv/delly/delly-0.6.1/*/Somatic_Germline_Quality_tagged.vcf',
        '/projects/POG/POG_data/(.*?)/sv/delly/delly-0.6.1/(.*?)/Somatic_Germline_Quality_tagged.vcf',
        ['patient', 'samples'],
    )
    delly_df = delly_df[delly_df.samples.str.split('_').apply(lambda x: len(x) == 2)] # filter out comparisons of >2 samples
    delly_df[['lib1', 'lib2']] = delly_df.samples.str.split('_', 1, expand=True)

    abyss_df = glob_to_df(
        '/projects/POG/POG_data/POG*/wgs/*/trans-ABySS/fusions/LSR.tsv',
        '/projects/POG/POG_data/(.*?)/wgs/(.*?)/trans-ABySS/fusions/LSR.tsv',
        ['patient', 'samples'],
    )
    abyss_df[['prefix', 'lib']] = abyss_df.samples.str.split('_t_', 1, expand=True)

    full_df = full_df[full_df.patient.isin(delly_df.patient) & full_df.patient.isin(abyss_df.patient)]
    full_df = full_df[full_df.t_lib.isin(delly_df.lib1) | full_df.t_lib.isin(delly_df.lib2)]
    full_df = full_df[full_df.t_lib.isin(abyss_df.lib)]
    return full_df

def resolve_delly_path(wildcards):
    lib_configurations = [
        '_'.join([wildcards.t_lib, wildcards.n_lib]),
        '_'.join([wildcards.n_lib, wildcards.t_lib])
    ]
    for configuration in lib_configurations:
        delly_versions = glob.glob('/projects/POG/POG_data/{patient}/sv/delly/delly-*'.format(patient = wildcards.patient))
        if len(delly_versions) == 0:
            raise(FileExistenceError('No DELLY file exists for POG patient {0} sample {1}_t_{2}_{3}_n_{4}'.format(
                wildcards.patient, wildcards.t_prefix, wildcards.t_lib, wildcards.n_prefix, wildcards.n_lib
            )))

        # There is some inconsistent DELLY file naming, so the next lines should generate every possible
        # file naming combination for these files.

        delly_filenames = [
            'Somatic_Germline_Quality_tagged.vcf',
            'Somatic_Germline_tagged.vcf',
            'Somatic_Quality_tagged.vcf',
            'Somatic_tagged.vcf'
        ]
        delly_filenames = \
            delly_filenames + \
            [d.lower() for d in delly_filenames]
        delly_filenames = delly_filenames + ['{}_{}'.format(wildcards.patient, d) for d in delly_filenames]

        for filename in delly_filenames:
            full_path = '{delly}/{libs}/{filename}'.format(
                delly = delly_versions[0],
                libs = configuration.strip(),
                filename = filename
            )
            if os.path.isfile(full_path):
                return([full_path])

    raise(FileExistenceError('No DELLY file exists for POG patient {0} sample {1}_t_{2}_{3}_n_{4}'.format(
        wildcards.patient, wildcards.t_prefix, wildcards.t_lib, wildcards.n_prefix, wildcards.n_lib
    )))

target_functions['pog'] = get_pog_targets

rule pog_snvs_and_indels:
    input:
        ancient(pog_project_root + '/{patient}/wgs/{sample}/strelka/bwa/results/passed.somatic.{mut_type}.vcf')
    output:
        'data/pog/patients/{patient}/{sample}/somatic_{mut_type}.vcf',
    wildcard_constraints:
        mut_type='(snvs|indels)' 
    shell:
        'cp -f {input} {output}'

rule pog_segments:
    input:
        ancient(pog_project_root + '/{patient}/wgs/{sample}/reviewed/loh/results/apolloh_out_segs.txt')
    output:
        'data/pog/patients/{patient}/{sample}/segments.tsv'
    log:
        'logs/data/pog/patients/{patient}/{sample}/segments.log'
    run:
        data = pd.read_csv(
            input[0],
            sep='\t',
            names=['chr', 'start', 'end', 'width', 'nvar', 'copy_number', 'lohtype', 'major_allele', 'minor_allele']
        )
        data.to_csv(output[0], sep='\t', index=False)

rule pog_sv:
    input:
        delly = resolve_delly_path,
        abyss = ancient('/projects/POG/POG_data/{patient}/wgs/{t_prefix}_t_{t_lib}/trans-ABySS/fusions/LSR.tsv'),
        match = ancient('/projects/POG/POG_data/{patient}/wgs/{t_prefix}_t_{t_lib}_{n_prefix}_n_{n_lib}')
    output:
        'data/pog/patients/{patient}/{t_prefix}_t_{t_lib}_{n_prefix}_n_{n_lib}/somatic_sv.tsv'
    log:
        'logs/pog/patients/{patient}/{t_prefix}_t_{t_lib}_{n_prefix}_n_{n_lib}/somatic_sv.log'
    shell:
        'Rscript scripts/projects/pog/merge_delly_abyss.R -d {input.delly} -a {input.abyss} -o {output}'

