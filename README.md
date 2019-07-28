# HRDetect

This implementation of the HRDetect pipeline is intended for RESEARCH USE ONLY and NOT FOR THE PURPOSE OF INFORMING CLINICAL DECISION-MAKING.

This repository contains an implementation of [HRDetect](https://www.nature.com/articles/nm.4292) used by Eric Y. Zhao at the Genome Sciences Centre.

HRDetect is run via a [Snakemake](http://snakemake.readthedocs.io/en/latest/) pipeline. It depends upon a working installation of R, with numerous dependencies. These may be installed by running `make dependencies`, which builds a miniconda environment with the necessary installations. To source that environment after it is installed, simply run `source dependencies/miniconda3/bin/activate dependencies` prior to running the pipeline.

Depending on your system and needs, some of these dependencies may require tweaking. Notably, the dependency installer assumes a linux-based environment.

Also necessary are two in-house tools called `SignIT` and `hrdtools`, which can be acquired by running `make`.

```
> make

if [ -d git/hrdtools ]; \
    then(cd git/hrdtools && git pull); \
    else git clone git@github.com:eyzhao/hrdtools.git git/hrdtools; \
    fi
Cloning into 'git/hrdtools'...
remote: Counting objects: 97, done.
remote: Compressing objects: 100% (6/6), done.
remote: Total 97 (delta 1), reused 0 (delta 0), pack-reused 91
Receiving objects: 100% (97/97), 163.64 KiB | 821.00 KiB/s, done.
Resolving deltas: 100% (16/16), done.
if [ -d git/SignIT ]; \
    then(cd git/SignIT && git pull); \
    else git clone git@github.com:eyzhao/SignIT.git git/SignIT; \
    fi
Cloning into 'git/SignIT'...
remote: Counting objects: 394, done.
remote: Compressing objects: 100% (113/113), done.
remote: Total 394 (delta 50), reused 58 (delta 22), pack-reused 259
Receiving objects: 100% (394/394), 1.11 MiB | 1.44 MiB/s, done.
Resolving deltas: 100% (153/153), done.
```

# Adding a project-specific pipeline

HRDetect expects data files of specific types and formats at specific locations. All files are housed at `data/{project}/{subproject}/patients/{patient}/{sample}/`, where `{patient}` and `{sample}` together form a unique ID pair for a given sample. Within each such directory, there must be four files specific to the sample.

1. `segments.tsv`
2. `somatic_indels.vcf`
3. `somatic_snvs.vcf`
4. `somatic_sv.tsv`

If you would like to use the Snakemake pipeline as is, then you can provide a project-specific file under the `projects` directory. Some projects files are already there as an example. You can then link to the project by adding a line `include: "project/myproject.smk"` in Snakefile.

If you would like to construct your own pipeline structure, please feel free to use the scripts in the `scripts` folder as needed.

## Example pipeline

For those interested in constructing their own working pipeline using Snakemake, there is an example project named `example` within this repository to demonstrate how this could work.

- An example of raw input data can be found in `examples/example`.
- This example data is copied and/or processed into the appropriate format, which is then stored in `data/example`. The Snakemake pipeline commands that pre-process the raw data are found in `projects/example.smk`.
- The expected output after running the full pipeline can be found in `output/example`
- The `Snakefile` is currently set up to run this example data. If you would like to run an integration test yourself, try deleting part or all of `output/example` and `data/example` and running `snakemake -p` (adding `-p` doesn't change how Snakemake runs, but will show you the commands being issued in each step).

Note that you may still run into errors/issues depending on your specific environment setup. As this pipeline has not yet been tested across a wide variety of environments, we are continuing to work out the kinks.

## segments.tsv

This is a file with segmented CNV/LOH calls with at least 5 columns.

1. `chr`: The chromosome name
2. `start`: Start position of CNV/LOH call
3. `end`: End position of CNV/LOH call
4. `copy_number`: The tumour copy number of the segment
5. `lohtype`: The type of LOH state. Should be amongst the following:
    - `ASCNA`: Allele-specific copy number amplification
    - `BCNA`: Balanced copy number amplification
    - `HET`: Heterozygous (normal)
    - `NLOH`: Neutral LOH (loss of heterozygosity, but 2 copies present)
    - `DLOH`: Deletion LOH
    - `ALOH`: Amplification LOH

## somatic_indels.vcf

A VCF file containing indels which can be parsed by R using the `readVCF()` function of [VariantAnnotation](http://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html).

## somatic_snvs.vcf

A VCF file containing SNVs which can be parsed by R using the `readVCF()` function of [VariantAnnotation](http://bioconductor.org/packages/release/bioc/html/VariantAnnotation.html).

## somatic_sv.tsv

A tab-delimited file with structural variant data. Should contain the following columns:

1. `chr1`: Name of the first chromosome involved
2. `pos1`: Coordinate of the SV breakpoint corresponding to the first chromosome
3. `chr2`: Name of the second chromosome involved
4. `pos2`: Coordinate of the SV breakpoint corresponding to the second chromosome
5. `type`: Can take on values `DEL`, `DUP`, `TRA`, or `INV`. Unless the value is `TRA`, the two chromosomes should be the same.

# Citation

If you use HRDtools or this implementation of HRDetect in your publication, please cite the following study:

[Zhao EY, Shen Y, Pleasance E, ... Jones SJM, Homologous Recombination Deficiency and Platinum-Based Therapy Outcomes in Advanced Breast Cancer. Clinical Cancer Research. 2017 Dec 15;23(24):7521-7530](https://www.ncbi.nlm.nih.gov/pubmed/29246904)


