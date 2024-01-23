#!/bin/bash

set -e

neurodocker generate singularity \
	--pkg-manager apt \
   	--base-image debian:bullseye-slim \
   	--install libblas-dev liblapack-dev git python3-dev python3-virtualenv python3-tk gcc g++ libeigen3-dev r-base libicu-dev libgmp-dev libmpfr-dev mpi-default-bin mpi-default-dev libcgal-dev gmsh libfreetype6-dev libxml2-dev libxslt-dev pip   \
   	--miniconda \
   		version=latest \
   		pip_install="Connectome-Spatial-Smoothing" \
		conda_install="python=3.10" \
	--run-bash 'git clone https://github.com/scilus/scilpy.git && mv scilpy /opt/ && pip install --upgrade pip && pip install -e /opt/scilpy/.' \
   	 > singularity_css

