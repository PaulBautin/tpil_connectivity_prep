#!/usr/bin/env nextflow
nextflow.enable.dsl=2

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["outlier_alpha":"$params.outlier_alpha",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info ""
log.info "TPIL Connectivity Prep Pipeline"
log.info "=================================="
log.info "Start time: $workflow.start"
log.info ""
log.info "[Input info]"
log.info "Input tractoflow folder: $params.input_tr"
log.info "Input freesurfer folder: $params.input_fs"
log.info "Input freesurfer licence file: $params.licence_fs"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}


process Subcortex_segmentation {
    input:
    tuple val(sid), path(T1nativepro_brain), path(affine), path(warp), path(inv_warp), file(t1_diffpro_brain)

    output:
    tuple val(sid), file("${sid}_all_fast_firstseg.nii.gz"), emit: sub_parcels
    tuple val(sid), file("${sid}*_first_diff.ply"), emit: sub_surfaces


    script:
    """
    run_first_all -i ${T1nativepro_brain} -o ${sid} -b -v
    #source activate env_scil
    #pip install --no-index 'nibabel<4'
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-L_Thal_first.vtk  ${sid}-L_Thal_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-L_Caud_first.vtk  ${sid}-L_Caud_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-L_Puta_first.vtk  ${sid}-L_Puta_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-L_Pall_first.vtk  ${sid}-L_Pall_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-BrStem_first.vtk  ${sid}-BrStem_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-L_Hipp_first.vtk  ${sid}-L_Hipp_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-L_Amyg_first.vtk  ${sid}-L_Amyg_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-L_Accu_first.vtk  ${sid}-L_Accu_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-R_Thal_first.vtk  ${sid}-R_Thal_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-R_Caud_first.vtk  ${sid}-R_Caud_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-R_Puta_first.vtk  ${sid}-R_Puta_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-R_Pall_first.vtk  ${sid}-R_Pall_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-R_Hipp_first.vtk  ${sid}-R_Hipp_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-R_Amyg_first.vtk  ${sid}-R_amyg_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-R_Accu_first.vtk  ${sid}-R_Accu_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    python /home/pabaua/dev_tpil/tpil_dmri/surface/first2diff.py ${sid}-R_Caud_first.vtk  ${sid}-R_Caud_first_diff.ply --ref_img ${T1nativepro_brain}  -f
    """
}

process Tractography_registration {
    input:
    tuple val(sid), path(T1nativepro_brain), path(affine), path(warp), path(tracto_pft)

    output:
    tuple val(sid), file("${bname}_t1.trk")

    script:
    bname = tracto_pft.name.split('.trk')[0]
    """
    source activate env_scil
    scil_apply_transform_to_tractogram.py ${tracto_pft} $T1nativepro_brain ${affine} ${bname}_t1.trk --in_deformation ${warp} -v --reverse_operation -f
    """
}


process Tractography_filtering {
    input:
    tuple val(sid), path(tracto_pft), path(seg_first)

    output:
    tuple val(sid), file("${bname}__NAc_proj.trk")

    script:
    bname = tracto_pft.name.split('.trk')[0]
    """
    source activate env_scil
    scil_filter_tractogram.py ${tracto_pft} ${bname}__NAc_proj.trk --atlas_roi ${seg_first} 26 any include -f -v
    """
}


process Fs_ciftify {
    input:
    tuple val(sid), path(seg_first), path(lh_surf), path(rh_surf), path(labels)

    output:
    tuple val(sid), file("${sid}.dscalar.nii"), emit: dscalar
    tuple val(sid), file("${sid}/*"), emit: dir
    tuple val(sid), file("${sid}/T1w/fsaverage_LR32k/${sid}.L.white.32k_fs_LR.surf.gii"), file("${sid}/T1w/fsaverage_LR32k/${sid}.R.white.32k_fs_LR.surf.gii"), emit: surfaces
    

    script:
    """   
    export PATH="/opt/fsl-6.0.5.1/bin:$PATH"
    export PATH="/opt/workbench/bin_linux64:$PATH"
    source activate env_cifti
    set +u
    source /opt/freesurfer-7.3.1/SetUpFreeSurfer.sh
    echo $PATH
    wb_command -volume-label-import  ${seg_first} ${labels} ${sid}_volume_label.nii.gz
    wb_command -cifti-create-dense-scalar ${sid}.dscalar.nii -volume ${seg_first} ${sid}_volume_label.nii.gz -left-metric $lh_surf -right-metric $rh_surf
    ciftify_recon_all --fs-subjects-dir $SUBJECTS_DIR --ciftify-work-dir . --resample-to-T1w32k $sid -v
    """
}


process Save_cifti {
    input:
    tuple val(sid),  path(tracto_pft), path(lh_surf), path(rh_surf), path(dscalar)

    output:
    tuple val(sid), file("${sid}.dscalar.nii"), emit: dscalar
    

    script:
    """   
    source activate env_scil
    python /home/pabaua/dev_tpil/tpil_dmri/surface/CSS_2.py $tracto_pft $lh_surf $rh_surf $dscalar ${sid}.dscalar.nii
    """
}


process BN_to_fs {
    memory_limit='6 GB'
    cpus=4

    input:
    tuple val(sid), path(T1nativepro_brain)

    output:
    tuple val(sid), file("lh.BN_Atlas.annot"), file("rh.BN_Atlas.annot"), file("BN_Atlas.nii.gz")

    script:
    """
    mris_ca_label -l $SUBJECTS_DIR/${sid}/label/lh.cortex.label ${sid} lh sphere.reg $SUBJECTS_DIR/lh.BN_Atlas.gcs lh.BN_Atlas.annot -t $SUBJECTS_DIR/BN_Atlas_210_LUT.txt
    cp $SUBJECTS_DIR/${sid}/label/lh.BN_Atlas.annot .
    mris_ca_label -l $SUBJECTS_DIR/${sid}/label/rh.cortex.label ${sid} rh sphere.reg $SUBJECTS_DIR/rh.BN_Atlas.gcs rh.BN_Atlas.annot -t $SUBJECTS_DIR/BN_Atlas_210_LUT.txt
    cp $SUBJECTS_DIR/${sid}/label/rh.BN_Atlas.annot .
    mri_aparc2aseg --s ${sid} --annot BN_Atlas --o $SUBJECTS_DIR/${sid}/mri/BN_Atlas.nii.gz
    cp $SUBJECTS_DIR/${sid}/mri/BN_Atlas.nii.gz .
    """
}


process Parcels_to_subject {
    memory_limit='6 GB'
    cpus=4

    input:
    tuple val(sid), path(fs_seg_BN), path(sub_seg)

    output:
    tuple val(sid), file("mask_subcortical.nii.gz"), file("mask_subcortical_dil.nii.gz"), file("fs_BN_dilated.nii.gz"), file("fs_BN_clipped.nii.gz"), file("BN_atlas_first.nii.gz"), file("BN_atlas_first_dil.nii.gz")

    script:
    """
    /opt/ants-2.3.2/bin/antsApplyTransforms -d 3 -i ${fs_seg_BN} -r ${sub_seg} -o fs_BN.nii.gz  -n genericLabel
    source activate env_scil
    scil_image_math.py convert fs_BN.nii.gz fs_BN.nii.gz --data_type uint16 -f
    scil_split_volume_by_labels.py fs_BN.nii.gz --scilpy_lut freesurfer_subcortical --out_dir subcortical_parcels
    scil_image_math.py union subcortical_parcels/* mask_subcortical.nii.gz -f
    scil_dilate_labels.py mask_subcortical.nii.gz mask_subcortical_dil.nii.gz  --distance 1.1 -f

    scil_image_math.py lower_threshold fs_BN.nii.gz 0.1 mask_brain.nii.gz -f
    scil_dilate_labels.py mask_brain.nii.gz mask_brain_dil.nii.gz  --distance 1.1 -f

    scil_image_math.py lower_clip fs_BN.nii.gz 2000 fs_BN_clipped_2000.nii.gz
    scil_image_math.py upper_clip fs_BN_clipped_2000.nii.gz 2211 fs_BN_clipped_2000.nii.gz -f
    scil_remove_labels.py fs_BN_clipped_2000.nii.gz  fs_BN_clipped_2000.nii.gz -i 2000 2211 -f
    scil_image_math.py subtraction fs_BN_clipped_2000.nii.gz 2000 fs_BN_clipped_2000.nii.gz --exclude_background --data_type int16 -f

    scil_image_math.py lower_clip fs_BN.nii.gz 1000 fs_BN_clipped_1000.nii.gz
    scil_image_math.py upper_clip fs_BN_clipped_1000.nii.gz 1211 fs_BN_clipped_1000.nii.gz -f
    scil_remove_labels.py fs_BN_clipped_1000.nii.gz  fs_BN_clipped_1000.nii.gz -i 1000 1211 -f
    scil_image_math.py subtraction fs_BN_clipped_1000.nii.gz 1000 fs_BN_clipped_1000.nii.gz --exclude_background --data_type int16 -f

    scil_combine_labels.py fs_BN_clipped.nii.gz --volume_ids fs_BN_clipped_2000.nii.gz all --volume_ids fs_BN_clipped_1000.nii.gz  all -f
    scil_dilate_labels.py fs_BN_clipped.nii.gz fs_BN_dilated.nii.gz  --distance 1.5 --mask mask_brain_dil.nii.gz

    #scil_image_math.py convert ${sub_seg} sub_seg.nii.gz --data_type uint16 -f
    #scil_image_math.py convert fs_BN_clipped.nii.gz fs_BN_clipped.nii.gz --data_type uint16 -f
    scil_image_math.py addition ${sub_seg} 1000 sub_seg_add_1000.nii.gz --exclude_background --data_type int16 -f
    scil_combine_labels.py  BN_atlas_first.nii.gz --volume_ids sub_seg_add_1000.nii.gz all --volume_ids fs_BN_clipped.nii.gz all -f
    scil_combine_labels.py  BN_atlas_first_dil.nii.gz --volume_ids sub_seg_add_1000.nii.gz all --volume_ids fs_BN_dilated.nii.gz all -f
    """
}

process Connectlow_prep {
    publishDir = {"./results/$sid"}

    input:
    tuple val(sid), path(tracto), path(labels_schaefer), path(labels_BN), path(T1diffpro_brain), path(dwi), path(bval), path(bvec), path(peaks), path(fodf)

    output:
    tuple val(sid), file("${sid}__tracking_pft.trk"), file("${sid}__schaefer_labels.nii.gz"), file("${sid}__BN_labels.nii.gz"), file("${sid}__t1.nii.gz"),file("${sid}__dwi.nii.gz"),file("${sid}__dwi.bval"),file("${sid}__dwi.bvec"), file("${peaks}"), file("${sid}__fodf.nii.gz")

    script:
    """
    cp ${tracto} ${sid}__tracking_pft.trk
    cp ${labels_schaefer} ${sid}__schaefer_labels.nii.gz
    cp ${labels_BN} ${sid}__BN_labels.nii.gz
    cp ${T1diffpro_brain} ${sid}__t1.nii.gz
    cp ${dwi} ${sid}__dwi.nii.gz
    cp ${bval} ${sid}__dwi.bval
    cp ${bvec} ${sid}__dwi.bvec
    """
}


process Accumbofrontal_seg {
    memory_limit='6 GB'

    input:
    tuple val(sid), path(tracto_L), path(tracto_R), path(labels_schaefer), path(labels_BN), path(brain_mask)

    output:
    tuple val(sid), file("${sid}__tracto_filtered_L.trk"), file("${sid}__tracto_cleaned_L.trk"), file("${sid}__tracto_filtered_R.trk"), file("${sid}__tracto_cleaned_R.trk")

    script:
    """
    source activate env_scil
    scil_image_math.py convert ${brain_mask} brain_mask.nii.gz --data_type uint16 -f
    scil_filter_tractogram.py ${tracto_L} ${sid}__tracto_filtered_L.trk --atlas_roi ${labels_BN} 27 any include -f --reference ${labels_BN}
    scil_cut_streamlines.py ${sid}__tracto_filtered_L.trk brain_mask.nii.gz ${sid}__tracto_filtered_L.trk --biggest_blob -f
    scil_outlier_rejection.py ${sid}__tracto_filtered_L.trk ${sid}__tracto_cleaned_L.trk --alpha 0.6 -f

    scil_filter_tractogram.py ${tracto_R} ${sid}__tracto_filtered_R.trk --atlas_roi ${labels_BN} 28 any include -f --reference ${labels_BN}
    scil_cut_streamlines.py ${sid}__tracto_filtered_R.trk brain_mask.nii.gz ${sid}__tracto_filtered_R.trk --biggest_blob -f
    scil_outlier_rejection.py ${sid}__tracto_filtered_R.trk ${sid}__tracto_cleaned_R.trk --alpha 0.6 -f
    """
}


process Accumbofrontal_seg_old {
    memory_limit='6 GB'

    input:
    tuple val(sid), path(tracto), path(labels_schaefer), path(labels_BN)

    output:
    tuple val(sid), file("${sid}__tracto_filtered_old_L_BN.trk"), file("${sid}__tracto_filtered_old_L_SC.trk"), emit:filtered
    tuple val(sid), file("${sid}__1026_47_L_cleaned_BN.trk"), emit:cleaned_BN
    tuple val(sid), file("${sid}__26_1111_L_cleaned_SC.trk"), emit:cleaned_SC

    script:
    """
    source activate env_scil
    scil_filter_tractogram.py ${tracto} ${sid}__tracto_filtered_old_L_BN.trk --atlas_roi ${labels_BN} 47 either_end include -f
    scil_filter_tractogram.py ${sid}__tracto_filtered_old_L_BN.trk ${sid}__tracto_filtered_old_L_BN.trk --atlas_roi ${labels_BN} 1026 either_end include -f
    scil_outlier_rejection.py ${sid}__tracto_filtered_old_L_BN.trk ${sid}__1026_47_L_cleaned_BN.trk --alpha 0.6 -f

    scil_filter_tractogram.py ${tracto} ${sid}__tracto_filtered_old_L_SC.trk --atlas_roi ${labels_schaefer} 1111 either_end include -f
    scil_filter_tractogram.py ${sid}__tracto_filtered_old_L_SC.trk ${sid}__tracto_filtered_old_L_SC.trk --atlas_roi ${labels_schaefer} 26 either_end include -f
    scil_outlier_rejection.py ${sid}__tracto_filtered_old_L_SC.trk ${sid}__26_1111_L_cleaned_SC.trk --alpha 0.6 -f
    """
}


process Register_Bundle {
    input:
    tuple val(sid), file(bundle), file(native_anat), file(template)

    output:
    tuple val(sid), file("${bname}_mni.trk"), emit: bundle

    script:
    bname = bundle.name.split('.trk')[0]
    """
    /opt/ants-2.3.2/bin/antsRegistrationSyNQuick.sh -d 3 -f ${native_anat} -m ${template} -t s -o ${sid}__output
    source activate env_scil
    scil_apply_transform_to_tractogram.py $bundle $template ${sid}__output0GenericAffine.mat ${bname}_mni.trk --in_deformation ${sid}__output1Warp.nii.gz --reverse_operation
    """
}

process Bundle_Pairwise_Comparaison_Inter_Subject {
    publishDir = {"./results_bundle/$task.process/$b_name"}
    input:
    tuple val(b_name), file(bundles)

    output:
    tuple val(b_name), file("${b_name}.json")

    script:
    """
    source activate env_scil
    scil_evaluate_bundles_pairwise_agreement_measures.py $bundles ${b_name}.json
    """
}

process Bundle_Pairwise_Comparaison_Intra_Subject {
    publishDir = {"./results_bundle/$task.process/$sid"}
    input:
    tuple val(sid), val(b_names), file(bundles)

    output:
    tuple val(sid), file("${b_names}_Pairwise_Comparaison.json")

    script:
    """
    source activate env_scil
    scil_evaluate_bundles_pairwise_agreement_measures.py $bundles ${b_names}_Pairwise_Comparaison.json
    """
}


process Compute_Centroid {
    input:
    tuple val(sid), file(bundle)

    output:
    tuple val(sid), file("${sid}__NAC_mPFC_L_centroid.trk")

    script:
    """
    source activate env_scil
    scil_compute_centroid.py ${bundle} ${sid}__NAC_mPFC_L_centroid.trk --nb_points 20
    """
}

process Bundle_QC_screenshot {
    input:
    tuple val(sid), file(bundle), file(ref_image)

    output:
    tuple val(sid), file("${sid}__${bname}.png")

    script:
    bname = bundle.name.split("__")[1].split('_L_')[0]
    """
    source activate env_scil
    scil_visualize_bundles_mosaic.py ${ref_image} ${bundle} ${sid}__${bname}.png -f --light_screenshot --no_information --ttf /usr/share/fonts/truetype/dejavu/DejaVuSans.ttf
    """
}



workflow {
    // Input files to fetch
    input_tractoflow = file(params.input_tr)
    input_freesurfer = file(params.input_fs)
    template = Channel.fromPath("$params.template")

    fs_brain = Channel.fromPath("$input_freesurfer/**/brain.mgz").map{[it.parent.parent.name, it.parent.parent.parent]}

    t1_nativepro_brain = Channel.fromPath("$input_tractoflow/*/Crop_T1/*__t1_bet_cropped.nii.gz").map{[it.parent.parent.name, it]}
    t1_diffpro_brain = Channel.fromPath("$input_tractoflow/*/Register_T1/*__t1_warped.nii.gz").map{[it.parent.parent.name, it]}
    t1_to_diff_affine = Channel.fromPath("$input_tractoflow/*/Register_T1/*__output0GenericAffine.mat").map{[it.parent.parent.name, it]}
    t1_to_diff_warp = Channel.fromPath("$input_tractoflow/*/Register_T1/*__output1Warp.nii.gz").map{[it.parent.parent.name, it]}
    t1_to_diff_inv_warp = Channel.fromPath("$input_tractoflow/*/Register_T1/*__output1InverseWarp.nii.gz").map{[it.parent.parent.name, it]}
    tracto_diff_pft = Channel.fromPath("$input_tractoflow/*/PFT_Tracking/*__pft_tracking_prob_wm_seed_0.trk").map{[it.parent.parent.name, it]}
    tracto_diff_L_Accu = Channel.fromPath("$input_tractoflow/*/Local_Tracking/*L*_tractogram.tck").map{[it.parent.parent.name, it]}
    tracto_diff_R_Accu = Channel.fromPath("$input_tractoflow/*/Local_Tracking/*R*_tractogram.tck").map{[it.parent.parent.name, it]}
    dwi_diff_pft = Channel.fromPath("$input_tractoflow/*/Resample_DWI/*__dwi_resampled.nii.gz").map{[it.parent.parent.name, it]}
    bval_diff_eddy = Channel.fromPath("$input_tractoflow/*/Eddy_Topup/*__bval_eddy").map{[it.parent.parent.name, it]}
    bvec_diff_eddy = Channel.fromPath("$input_tractoflow/*/Eddy_Topup/*__dwi_eddy_corrected.bvec").map{[it.parent.parent.name, it]}
    peaks_diff = Channel.fromPath("$input_tractoflow/*/FODF_Metrics/*__peaks.nii.gz").map{[it.parent.parent.name, it]}
    fodf_diff = Channel.fromPath("$input_tractoflow/*/FODF_Metrics/*__fodf.nii.gz").map{[it.parent.parent.name, it]}

    // hcp files
    lh_surf = Channel.fromPath('/home/pabaua/dev_tpil/data/Freesurfer/22-09-21_t1_clbp_hcp_output/zz_templates/L.atlasroi.32k_fs_LR.shape.gii')
    rh_surf = Channel.fromPath('/home/pabaua/dev_tpil/data/Freesurfer/22-09-21_t1_clbp_hcp_output/zz_templates/R.atlasroi.32k_fs_LR.shape.gii')
    labels = Channel.fromPath('/home/pabaua/dev_hiball/css_test/label_names.txt')



    main:
    // Subcortex segmentation with first + registration to diffusion space
    t1_nativepro_brain.combine(t1_to_diff_affine, by:0).combine(t1_to_diff_warp, by:0).combine(t1_to_diff_inv_warp, by:0).combine(t1_diffpro_brain, by:0).set{data_sub_seg}
    Subcortex_segmentation(data_sub_seg)

    t1_nativepro_brain.combine(t1_to_diff_affine, by:0).combine(t1_to_diff_warp, by:0).combine(tracto_diff_pft, by:0).set{data_tracto_reg}
    Tractography_registration(data_tracto_reg)

    Tractography_registration.out.combine(Subcortex_segmentation.out.sub_parcels, by:0).set{data_tracto_filt}
    Tractography_filtering(data_tracto_filt)

    Subcortex_segmentation.out.sub_parcels.combine(lh_surf).combine(rh_surf).combine(labels).set{data_ciftify}
    Fs_ciftify(data_ciftify)

    Tractography_registration.out.combine(Fs_ciftify.out.surfaces, by:0).combine(Fs_ciftify.out.dscalar, by:0).set{data_save_cifti}
    Save_cifti(data_save_cifti)

    // Register Atlas parcels into freesurfer space
    //BN_to_fs(t1_nativepro_brain)

    //Parcels_to_subject(BN_to_fs.out.map{[it[0], it[3]]}.combine(Subcortex_segmentation.out.sub_parcels, by:0))
}

