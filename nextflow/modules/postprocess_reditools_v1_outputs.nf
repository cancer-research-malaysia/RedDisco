// Reditools v1 Post-processing 
process POSTPROCESS_REDITOOLS_V1_OUTPUTS {
    tag "$sampleId"
    label 'postprocessReditools'
    container "${params.container__reditools}"
    publishDir "${params.outputDir}/reditools-v1/${sampleId}/postprocessed", mode: 'copy'

    input:
    // Receives the existing tuple structure
    tuple val(sampleId), path(bamFile), path(indexFile), path(dnaBamFile), path(dnaIndexFile)
    tuple val(_sample), path(reditoolsOutputDir)
    path genomeFa
    path spliceSitesAnnot
    path excludedContigs
    path rmskGtf
    path rmskGtfIndex
    path snpGtf
    path snpGtfIndex
    path rediportalsDbGtf
    path rediportalsDbGtfIndex

    output:
    // Pass the intermediate files needed for the final wrangle
    tuple val(sampleId), path("${sampleId}-knownEditing"), path("${sampleId}--pos.txt"), path("${sampleId}--pos-ALU.txt"), emit: for_wrangling
    tuple val(sampleId), path("${sampleId}-firstalu"), path("${sampleId}-second"), emit: out_tables
    // Additionally, publish all post-processed files
    path "outTable_*"
    path "${sampleId}-first"
    
    script:
    """
    # print DNA-RNA editing output directory for verification
    echo "[Reditools Postprocessing] DNA-RNA editing output directory: ${reditoolsOutputDir}"
    echo "[Reditools Postprocessing] Starting post-processing for raw RNA editing calls for sample ${sampleId}."
    # now run the postprocessing script
    if bash reditools-output-postprocessing.sh ${bamFile} ${reditoolsOutputDir} ${genomeFa} ${spliceSitesAnnot} ${excludedContigs} ${rmskGtf} ${snpGtf} ${rediportalsDbGtf} ${task.cpus} ${sampleId}; then
       echo "[Reditools Postprocessing] Post-processing completed successfully for sample ${sampleId}."
    else
       echo "[Reditools Postprocessing] Error: Post-processing failed for sample ${sampleId}." >&2
       exit 1
    fi


    """
}
