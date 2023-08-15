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
log.info "Input MNI template: $params.template"
log.info "Input freesurfer licence file: $params.licence_fs"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}


process Subcortex_segmentation {
    input:
    tuple val(sid), path(T1nativepro_brain), path(affine), path(warp), file(t1_diffpro_brain)

    output:
    tuple val(sid), file("${sid}__first_atlas_transformed.nii.gz"), emit: sub_parcels
    tuple val(sid), file("*.vtk"), file("*.nii.gz"), emit: sub_surfaces

    script:
    """
    run_first_all -i ${T1nativepro_brain} -o ${sid} -b -v
    /opt/ants-2.3.2/bin/antsApplyTransforms -d 3 -i ${sid}_all_fast_firstseg.nii.gz -t ${warp} -t ${affine} -r ${t1_diffpro_brain} -o ${sid}__first_atlas_transformed.nii.gz -n genericLabel
    """
}


process Atlas_to_fs {
    memory_limit='6 GB'
    cpus=4

    input:
    tuple val(sid), path(T1nativepro_brain)

    output:
    tuple val(sid), file("lh.BN_Atlas.annot"), file("lh.BN_Atlas.nii.gz"), file("rh.BN_Atlas.annot"), file("rh.BN_Atlas.nii.gz")

    script:
    """
    tkregister2_cmdl --mov $SUBJECTS_DIR/${sid}/mri/brain.mgz --noedit --s ${sid} --regheader --reg register.dat
    mris_ca_label -l $SUBJECTS_DIR/${sid}/label/lh.cortex.label ${sid} lh sphere.reg $SUBJECTS_DIR/lh.BN_Atlas.gcs lh.BN_Atlas.annot -t $SUBJECTS_DIR/BN_Atlas_210_LUT.txt
    cp $SUBJECTS_DIR/${sid}/label/lh.BN_Atlas.annot .
    mris_ca_label -l $SUBJECTS_DIR/${sid}/label/rh.cortex.label ${sid} rh sphere.reg $SUBJECTS_DIR/rh.BN_Atlas.gcs rh.BN_Atlas.annot -t $SUBJECTS_DIR/BN_Atlas_210_LUT.txt
    cp $SUBJECTS_DIR/${sid}/label/rh.BN_Atlas.annot .
    mri_label2vol --subject ${sid} --hemi lh --annot BN_Atlas --o $SUBJECTS_DIR/${sid}/mri/lh.BN_Atlas.nii.gz --temp $SUBJECTS_DIR/${sid}/mri/brain.mgz --reg register.dat --proj frac 0 1 0.01
    cp $SUBJECTS_DIR/${sid}/mri/lh.BN_Atlas.nii.gz .
    mri_label2vol --subject ${sid} --hemi rh --annot BN_Atlas --o $SUBJECTS_DIR/${sid}/mri/rh.BN_Atlas.nii.gz --temp $SUBJECTS_DIR/${sid}/mri/brain.mgz --reg register.dat --proj frac 0 1 0.01
    cp $SUBJECTS_DIR/${sid}/mri/rh.BN_Atlas.nii.gz .
    """
}


process Parcels_to_subject {
    memory_limit='6 GB'
    cpus=4

    input:
    tuple val(sid), path(fs_seg_lh), path(fs_seg_rh), path(sub_seg), file(t1_diffpro_brain)

    output:
    tuple val(sid), file("${sid}__fsatlas_transformed.nii.gz"), emit: fs_parcels
    tuple val(sid), file("${sid}__nativepro_seg_all.nii.gz"), emit: seg_all

    script:
    """
    #mri_convert $SUBJECTS_DIR/${sid}/mri/brainmask.mgz mask_brain.nii.gz
    mri_convert $SUBJECTS_DIR/${sid}/mri/ribbon.mgz mask_brain.nii.gz
    source activate env_scil
    scil_image_math.py lower_threshold mask_brain.nii.gz 1 mask_brain_bin.nii.gz
    scil_combine_labels.py out_labels.nii.gz --volume_ids ${fs_seg_lh} all --volume_ids ${fs_seg_rh} all
    scil_dilate_labels.py out_labels.nii.gz fs_labels_dilated.nii.gz --distance 1.5 --mask mask_brain_bin.nii.gz
    /opt/ants-2.3.2/bin/antsRegistrationSyN -d 3 -f ${t1_diffpro_brain} -m $SUBJECTS_DIR/${sid}/mri/brain.mgz -t s -o ${sid}__output
    /opt/ants-2.3.2/bin/antsApplyTransforms -d 3 -i fs_labels_dilated.nii.gz -t ${sid}__output1Warp.nii.gz -t ${sid}__output0GenericAffine.mat -r ${t1_diffpro_brain} -o ${sid}__fsatlas_transformed.nii.gz -n GenericLabel
    scil_image_math.py addition ${sub_seg} 1000 sub_seg_add_1000.nii.gz --exclude_background --data_type int16
    scil_image_math.py addition ${sid}__fsatlas_transformed.nii.gz 0 cortex_seg_add_0.nii.gz --exclude_background --data_type int16
    scil_combine_labels.py ${sid}__nativepro_seg_all.nii.gz --volume_ids cortex_seg_add_0.nii.gz all --volume_ids sub_seg_add_1000.nii.gz all
    """
}


process Connectlow_prep {
    publishDir = {"./results/$sid"}

    input:
    tuple val(sid), path(tracto), path(labels), path(T1diffpro_brain), path(dwi), path(bval), path(bvec), path(peaks), path(fodf)

    output:
    tuple val(sid), file("${sid}__tracking_pft.trk"), file("${sid}__labels.nii.gz"), file("${sid}__t1.nii.gz"),file("${sid}__dwi.nii.gz"),file("${sid}__dwi.bval"),file("${sid}__dwi.bvec"), file("${peaks}"), file("${sid}__fodf.nii.gz")

    script:
    """
    cp ${tracto} ${sid}__tracking_pft.trk
    cp ${labels} ${sid}__labels.nii.gz
    cp ${T1diffpro_brain} ${sid}__t1.nii.gz
    cp ${dwi} ${sid}__dwi.nii.gz
    cp ${bval} ${sid}__dwi.bval
    cp ${bvec} ${sid}__dwi.bvec
    """
}



workflow {
    // Input files to fetch
    input_tractoflow = file(params.input_tr)
    input_freesurfer = file(params.input_fs)

    fs_brain = Channel.fromPath("$input_freesurfer/**/brain.mgz").map{[it.parent.parent.name, it.parent.parent.parent]}

    t1_nativepro_brain = Channel.fromPath("$input_tractoflow/*/Crop_T1/*__t1_bet_cropped.nii.gz").map{[it.parent.parent.name, it]}
    t1_diffpro_brain = Channel.fromPath("$input_tractoflow/*/Register_T1/*__t1_warped.nii.gz").map{[it.parent.parent.name, it]}
    t1_to_diff_affine = Channel.fromPath("$input_tractoflow/*/Register_T1/*__output0GenericAffine.mat").map{[it.parent.parent.name, it]}
    t1_to_diff_warp = Channel.fromPath("$input_tractoflow/*/Register_T1/*__output1Warp.nii.gz").map{[it.parent.parent.name, it]}
    tracto_diff_pft = Channel.fromPath("$input_tractoflow/*/PFT_Tracking/*__pft_tracking_prob_wm_seed_0.trk").map{[it.parent.parent.name, it]}
    dwi_diff_pft = Channel.fromPath("$input_tractoflow/*/Resample_DWI/*__dwi_resampled.nii.gz").map{[it.parent.parent.name, it]}
    bval_diff_eddy = Channel.fromPath("$input_tractoflow/*/Eddy_Topup/*__bval_eddy").map{[it.parent.parent.name, it]}
    bvec_diff_eddy = Channel.fromPath("$input_tractoflow/*/Eddy_Topup/*__dwi_eddy_corrected.bvec").map{[it.parent.parent.name, it]}
    peaks_diff = Channel.fromPath("$input_tractoflow/*/FODF_Metrics/*__peaks.nii.gz").map{[it.parent.parent.name, it]}
    fodf_diff = Channel.fromPath("$input_tractoflow/*/FODF_Metrics/*__fodf.nii.gz").map{[it.parent.parent.name, it]}



    main:
    // Subcortex segmentation with first + registration to diffusion space
    t1_nativepro_brain.combine(t1_to_diff_affine, by:0).combine(t1_to_diff_warp, by:0).combine(t1_diffpro_brain, by:0).set{data_sub_seg}
    Subcortex_segmentation(data_sub_seg)

    // Register Atlas parcels into freesurfer space
    Atlas_to_fs(t1_nativepro_brain)

    // Apply and combine cortex and sub-cortex parcels
    Atlas_to_fs.out.map{[it[0], it[2], it[4]]}.combine(Subcortex_segmentation.out.sub_parcels, by:0).combine(t1_diffpro_brain, by:0).set{data_atlas_to_fs}
    Parcels_to_subject(data_atlas_to_fs)

    // Connectlow prep
    tracto_diff_pft.combine(Parcels_to_subject.out.seg_all, by:0).combine(t1_diffpro_brain, by:0).combine(dwi_diff_pft, by:0).combine(bval_diff_eddy, by:0).combine(bvec_diff_eddy, by:0).combine(peaks_diff, by:0).combine(fodf_diff, by:0).set{data_connectflow_prep}
    Connectlow_prep(data_connectflow_prep)

}

