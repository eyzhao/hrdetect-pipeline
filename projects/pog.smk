pog_project_root = '/projects/POG/POG_data'

class FileExistenceError(Exception):
    pass

def is_integer(string):
    try:
        int(string)
        return True
    except ValueError:
        return False

def collapse_strelka_runs(wildcards):
    ''' Given a POG sample with multiple strelka runs, this function automatically
        selects the most recent strelka run (having higher run index) and returns
        only that one for use as snakemake input.
    '''
    path_template = pog_project_root + '/{patient}/wgs/{t_prefix}_t_{t_lib}_{n_prefix}_n_{n_lib}/{t_lib}_{n_lib}/strelka/{{strelka_run_id}}/bwa/results/passed.somatic.{mut_type}.vcf'
    path = path_template.format(
        patient = wildcards.patient,
        t_prefix = wildcards.t_prefix,
        t_lib = wildcards.t_lib,
        n_prefix = wildcards.n_prefix,
        n_lib = wildcards.n_lib,
        mut_type = wildcards.mut_type
    )

    (all_strelka_runs,) = glob_wildcards(path)
    all_strelka_runs = [i for i in all_strelka_runs if is_integer(i)]
    numeric_strelka_runs = [int(i) for i in all_strelka_runs]
    chosen_run = [str(max(numeric_strelka_runs))]

    return expand(path, strelka_run_id = chosen_run)

def get_pog_snv_targets(wildcards):
    return glob_to_df(
        '/projects/POG/POG_data/POG???/wgs/*/*/strelka/*/bwa/results/passed.somatic.snvs.vcf',
        '/projects/POG/POG_data/(.*?)/wgs/(.*?)/.*?/strelka/.*?/bwa/results/passed.somatic.snvs.vcf',
        ['patient', 'samples'],
    )

def get_pog_microhomology_targets(wildcards):
    return glob_to_df(
        '/projects/POG/POG_data/POG???/wgs/*/*/strelka/*/bwa/results/passed.somatic.indels.vcf',
        '/projects/POG/POG_data/(.*?)/wgs/(.*?)/.*?/strelka/.*?/bwa/results/passed.somatic.indels.vcf',
        ['patient', 'samples'],
    )

def get_pog_hrd_score_targets(wildcards):
    return glob_to_df(
        '/projects/POG/POG_data/POG???/wgs/*/reviewed/loh/results/apolloh_out_segs.txt',
        '/projects/POG/POG_data/(.*?)/wgs/(.*?)/reviewed/loh/results/apolloh_out_segs.txt',
        ['patient', 'samples'],
    )

def get_pog_sv_targets(wildcards):
    full_df = glob_to_df(
        '/projects/POG/POG_data/POG???/wgs/*_t_*_*_n_*',
        '/projects/POG/POG_data/(.*?)/wgs/(.*)',
        ['patient', 'samples'],
    )
    full_df[['t_prefix', 't', 't_lib', 'n_prefix', 'n', 'n_lib']] = full_df.samples.str.split('_', expand=True)

    delly_df = glob_to_df(
        '/projects/POG/POG_data/POG???/sv/delly/delly-0.6.1/*/Somatic_Germline_Quality_tagged.vcf',
        '/projects/POG/POG_data/(.*?)/sv/delly/delly-0.6.1/(.*?)/Somatic_Germline_Quality_tagged.vcf',
        ['patient', 'samples'],
    )
    delly_df = delly_df[delly_df.samples.str.split('_').apply(lambda x: len(x) == 2)] # filter out comparisons of >2 samples
    delly_df[['lib1', 'lib2']] = delly_df.samples.str.split('_', 1, expand=True)

    abyss_df = glob_to_df(
        '/projects/POG/POG_data/POG???/wgs/*/trans-ABySS/fusions/LSR.tsv',
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
        full_path = '/projects/POG/POG_data/{patient}/sv/delly/delly-0.6.1/{libs}/Somatic_Germline_Quality_tagged.vcf'.format(
            patient = wildcards.patient,
            libs = configuration.strip()
        )
        if os.path.isfile(full_path):
            return([full_path])

    raise(FileExistenceError('No DELLY file exists for POG patient {0} sample {1}_t_{2}_{3}_n_{4}'.format(
        wildcards.patient, wildcards.t_prefix, wildcards.t_lib, wildcards.n_prefix, wildcards.n_lib
    )))

target_functions['pog'] = {
    'all': {
        'snv_signatures': get_pog_snv_targets,
        'sv_signatures': get_pog_sv_targets,
        'hrd_scores': get_pog_hrd_score_targets,
        'microhomology_scores': get_pog_microhomology_targets,
    }
}

rule pog_snvs_and_indels:
    input:
        ancient(collapse_strelka_runs)
    output:
        'data/pog/all/patients/{patient}/{t_prefix}_t_{t_lib}_{n_prefix}_n_{n_lib}/somatic_{mut_type}.vcf',
    log:
        'logs/pog/all/patients/{patient}/{t_prefix}_t_{t_lib}_{n_prefix}_n_{n_lib}/somatic_{mut_type}.log',
    wildcard_constraints:
        mut_type='(snvs|indels)' 
    shell:
        'cp -f {input} {output}'

rule pog_segments:
    input:
        ancient(pog_project_root + '/{patient}/wgs/{sample}/reviewed/loh/results/apolloh_out_segs.txt')
    output:
        'data/pog/all/patients/{patient}/{sample}/segments.tsv'
    log:
        'logs/data/pog/all/patients/{patient}/{sample}/segments.log'
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
        'data/pog/all/patients/{patient}/{t_prefix}_t_{t_lib}_{n_prefix}_n_{n_lib}/somatic_sv.tsv'
    log:
        'logs/pog/all/patients/{patient}/{t_prefix}_t_{t_lib}_{n_prefix}_n_{n_lib}/somatic_sv.log'
    shell:
        'Rscript ' + PROJECT_DIR + '/scripts/projects/pog/merge_delly_abyss.R -d {input.delly} -a {input.abyss} -o {output}'

