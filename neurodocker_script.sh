#!/bin/bash

set -e

neurodocker generate singularity \
	--pkg-manager apt \
   	--base-image ubuntu:focal-20230624 \
   	--arg DEBIAN_FRONTEND=noninteractive \
   	--arg DEBCONF_NONINTERACTIVE_SEEN=true \
   	--arg TZ=America/Toronto \
   	--install tzdata \
   	--user nonroot \
   	--install libblas-dev liblapack-dev git python3-dev python3-virtualenv python3-tk gcc g++ libeigen3-dev r-base libicu-dev libgmp-dev libmpfr-dev mpi-default-bin mpi-default-dev libcgal-dev gmsh libfreetype6-dev libxml2-dev libxslt-dev pip  \
   	--miniconda \
   		version=latest \
   		env_name='env_scil' \
   		env_exists=false \
   		conda_install="python=3.10" \
   	--run-bash 'source activate env_scil && git clone https://github.com/scilus/scilpy.git && mv scilpy /opt/ && pip install --upgrade pip && pip install -e /opt/scilpy/. && source deactivate' \
   	--ants version=2.3.2 \
	--freesurfer version=7.3.1 \
   	> singularity_prep

