dependencies_repo_url := git@github.com:eyzhao/bio-pipeline-dependencies.git
project_root := /projects/ezhao_prj/analyses/hrdetect-pipeline

###########################
### Directory Structure ###
###########################

meta:
	mkdir -p meta

paths:
	mkdir -p paths

#################
### Load Code ###
#################

dependencies:
	if [ -d $@ ]; \
	then (cd $@ && git checkout hrdetect && git pull); \
	else git clone ${dependencies_repo_url} $@ && cd $@ && git checkout hrdetect && git pull; \
	fi && \
	make

scripts/SignIT:
	if [ -d $@ ]; \
	then(cd $@ && git pull); \
	else git clone git@github.com:eyzhao/SignIT.git $@; \
	fi

###################################
### Load paths to all POG files ###
###################################

paths/flatfiles.txt: paths
	ls /projects/POG/POG_data/POG???/flatfile/POG???.tab > paths/flatfiles.txt

meta/all_pog_metadata.txt: paths/flatfiles.txt scripts/pipelines meta
	Rscript scripts/pipelines/pog-files/aggregate_metadata.R -p $< -o $@

paths/delly_paths.txt: paths
	ls /projects/POG/POG_data/POG*/sv/delly/*/*/Somatic_Germline_Quality_tagged.vcf > $@

paths/abyss_paths.txt: paths
	ls /projects/POG/POG_data/POG*/wgs/*/trans-ABySS/fusions/LSR.tsv > $@

paths/delly_abyss_paths.txt: paths/delly_paths.txt paths/abyss_paths.txt meta/all_pog_metadata.txt
	Rscript scripts/pipelines/sv-signatures/pair_sv_abyss_paths.R -s $< -a paths/abyss_paths.txt -m meta/all_pog_metadata.txt -o $@

paths/pog_merged_paths.txt: paths
	Rscript scripts/pipelines/pog-files/get_paths.R  \
		-G '/projects/POG/POG_data/POG*/flatfile/*.tab'  \
		-r '.*?\/(POG\d+)\.tab'  \
		-n 'pog_id'  \
		-p 'flatfile'  \
		-m 'pog_id' \
	| Rscript scripts/pipelines/pog-files/get_paths.R  \
		-G '/projects/POG/POG_data/POG*/wgs/*/reviewed/loh/results/apolloh_out_segs.txt'  \
		-r '.*?POG_data\/(POG\d+)\/wgs\/(.*?)\/.*'  \
		-n 'pog_id,comparison'  \
		-p 'loh_segs'  \
		-m 'pog_id' \
	| Rscript scripts/pipelines/pog-files/get_paths.R  \
		-G '/projects/POG/POG_data/POG*/wgs/*/reviewed/cnv/w200/w200_segs.txt'  \
		-r '.*?POG_data\/(POG\d+)\/wgs\/(.*?)\/.*'  \
		-n 'pog_id,comparison'  \
		-p 'cnv_segs'  \
		-m 'pog_id,comparison' \
	| Rscript scripts/pipelines/pog-files/get_paths.R  \
		-G '/projects/POG/POG_data/*/wgs/*/*/strelka/*/bwa/results/passed.somatic.snvs.vcf'  \
		-r '.*POG_data\/(POG\d+)\/wgs\/(.*?)\/.*?strelka\/(.*?)\/bwa.*'  \
		-n 'pog_id,comparison,strelka_run_id'  \
		-p 'snv'  \
		-m 'pog_id,comparison' \
	| Rscript -e "library(tidyverse); read_tsv(file('stdin')) %>% gather(file, path, -pog_id, -comparison, -strelka_run_id) %>% write_tsv('$@')"

