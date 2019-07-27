dependencies_repo_url := https://github.com/eyzhao/bio-pipeline-dependencies.git
project_root := /projects/ezhao_prj/analyses/hrdetect-pipeline

###########################
### Directory Structure ###
###########################

all: \
	git/hrdtools \
	git/SignIT

#################
### Load Code ###
#################

dependencies:
	if [ -d $@ ]; \
	then (cd $@ && git checkout hrdetect && git pull); \
	else git clone ${dependencies_repo_url} $@ && cd $@ && git checkout hrdetect && git pull; \
	fi && \
	make

git/SignIT:
	if [ -d $@ ]; \
	then(cd $@ && git pull); \
	else git clone https://github.com/eyzhao/SignIT.git $@; \
	fi

git/hrdtools:
	if [ -d $@ ]; \
	then(cd $@ && git pull); \
	else git clone https://github.com/eyzhao/hrdtools.git $@; \
	fi
