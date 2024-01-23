#!/bin/bash

# This would run the TPIL Connectictivity prep Pipeline with the following resources:


# Container image based on neurodocker
my_singularity_img='/home/pabaua/dev_scil/containers/fmriprep-23.2.0.simg' # or .sif
# main file for tpil_connectivity_prep
my_main_nf='/home/pabaua/dev_tpil/tpil_connectivity_prep/main_accumbofrontal_3.nf'
# Results folder of tractoflow
my_input_tr='/home/pabaua/dev_tpil/results/results_tracto/23-09-01_tractoflow_bundling'
# Results folder of freesurfer
my_input_fs='/home/pabaua/dev_tpil/data/Freesurfer/22-09-21_t1_clbp_freesurfer_output'
# Results folder of fmriprep
my_input_fp='/home/pabaua/dev_tpil/results/results_fmriprep/24-01-16_fmriprep'
# MNI template masked (the same as commonly used by the scil)
my_template='/home/pabaua/dev_scil/atlas/mni_masked.nii.gz'
# Freesurfer licence file
my_licence_fs='/home/pabaua/dev_tpil/data/Freesurfer/license.txt'


nextflow run $my_main_nf  \
  --input_tr $my_input_tr \
  --input_fs $my_input_fs \
  --input_fp $my_input_fp \
  --template $my_template \
  --licence_fs $my_licence_fs \
  --template $my_template \
  -with-singularity $my_singularity_img \
  -profile local \
  -resume \
  --cleanenv
