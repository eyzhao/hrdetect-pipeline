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
