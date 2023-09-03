#!/bin/bash

# This would run the TPIL Bundle Segmentation Pipeline with the following resources:
#     Prebuild Singularity images: https://scil.usherbrooke.ca/pages/containers/
#     Brainnetome atlas in MNI space: https://atlas.brainnetome.org/download.html
#     FA template in MNI space: https://brain.labsolver.org/hcp_template.html


#SBATCH --nodes=1              # --> Generally depends on your nb of subjects.
                               # See the comment for the cpus-per-task. One general rule could be
                               # that if you have more subjects than cores/cpus (ex, if you process 38
                               # subjects on 32 cpus-per-task), you could ask for one more node.
#SBATCH --cpus-per-task=32     # --> You can see here the choices. For beluga, you can choose 32, 40 or 64.
                               # https://docs.computecanada.ca/wiki/B%C3%A9luga/en#Node_Characteristics
#SBATCH --mem=0                # --> 0 means you take all the memory of the node. If you think you will need
                               # all the node, you can keep 0.
#SBATCH --time=48:00:00

#SBATCH --mail-user=paul.bautin@polymtl.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL


module load StdEnv/2020 java/14.0.2 nextflow/22.10.8 apptainer/1.1.8


my_singularity_img='/home/pabaua/projects/def-pascalt-ab/pabaua/dev_scil/containers/singularity_container.sif' # or .sif
my_main_nf='/home/pabaua/projects/def-pascalt-ab/pabaua/dev_tpil/tpil_connectivity_prep/main_accumbofrontal.nf'
my_input_tr_con='/home/pabaua/scratch/tpil_dev/results/control/23-07-05_tractoflow_bundling/results'
my_input_fs_con='/home/pabaua/scratch/tpil_dev/data/23_02_09_control_freesurfer_output'
my_input_tr_clbp='/home/pabaua/scratch/tpil_dev/results/clbp/23-07-05_tractoflow_bundling/results'
my_input_fs_clbp='/home/pabaua/scratch/tpil_dev/data/22-09-21_t1_clbp_freesurfer_output'
my_licence_fs='/home/pabaua/projects/def-pascalt-ab/pabaua/dev_tpil/tpil_connectivity_prep/freesurfer_data/license.txt'
my_template='/home/pabaua/projects/def-pascalt-ab/pabaua/dev_scil/atlas/mni_masked.nii.gz'


nextflow run $my_main_nf  \
  --input_tr $my_input_tr_con \
  --input_fs $my_input_fs_con \
  --licence_fs $my_licence_fs \
  --template $my_template \
  -with-singularity $my_singularity_img \
  -profile compute_canada \
  -resume


# module load StdEnv/2020 java/14.0.2 nextflow/21.10.3 apptainer/1.1.8


# my_singularity_img='/home/pabaua/projects/def-pascalt-ab/pabaua/dev_scil/containers/scilus_1.5.0.sif' # or .img
# my_main_nf='/home/pabaua/projects/def-pascalt-ab/pabaua/dev_scil/connectoflow/main.nf'
# my_input_con_BN='/home/pabaua/scratch/tpil_dev/results/control/23-08-17_connectflow/results'
# my_input_clbp_BN='/home/pabaua/scratch/tpil_dev/results/clbp/23-08-17_connectflow/results'
# my_input_con_schaefer='/home/pabaua/scratch/tpil_dev/results/control/23-08-17_connectflow_schaefer/results'
# my_input_clbp_schaefer='/home/pabaua/scratch/tpil_dev/results/clbp/23-08-17_connectflow_schaefer/results'
# my_template='/home/pabaua/projects/def-pascalt-ab/pabaua/dev_scil/atlas/mni_masked.nii.gz'
# my_labels_list_BN='/home/pabaua/projects/def-pascalt-ab/pabaua/dev_tpil/tpil_connectivity_prep/freesurfer_data/atlas_brainnetome_first_label_list.txt'
# my_labels_list_schaefer='/home/pabaua/projects/def-pascalt-ab/pabaua/dev_tpil/tpil_connectivity_prep/freesurfer_data/atlas_schaefer_200_first_label_list.txt'

# NXF_DEFAULT_DSL=1 nextflow run $my_main_nf \
#   --input $my_input_con_BN \
#   --labels_list $my_labels_list_BN \
#   --labels_img_prefix 'BN_' \
#   --template $my_template \
#   --apply_t1_labels_transfo false \
#   -with-singularity $my_singularity_img \
#   -resume

