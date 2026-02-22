// Reditools v1 Post-processing (rewritten) 
process POSTPROCESS_REDITOOLS_V1_OUTPUTS {
    tag "$sampleId"
    label 'postprocessReditools'
    container "${params.container__reditools}"
    publishDir "${params.outputDir}/reditools-v1/${sampleId}/postprocessed", mode: 'copy'

    input:
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
    // Direct equivalent of the old for_wrangling emit:
    //   ${SAMPLE_ID}--pos.txt     -> ${sampleId}-NONALU-NONREP--novelEditing.txt
    //   ${SAMPLE_ID}--pos-ALU.txt -> ${sampleId}-ALU--novelEditing.txt
    tuple val(sampleId),
          path("${sampleId}-knownEditing"),
          path("${sampleId}-NONALU-NONREP--novelEditing.txt"),
          path("${sampleId}-ALU--novelEditing.txt"),
          emit: for_wrangling

    // REDItools re-analysis output directories (novel ALU and novel NONALU+NONREP)
    tuple val(sampleId),
          path("${sampleId}-firstalu"),
          path("${sampleId}-second"),
          emit: out_tables

    // Final labeled outputs and summary statistics
    tuple val(sampleId),
          path("final_output_files/${sampleId}-allEditing.tsv"),
          path("final_output_files/${sampleId}-AG-Subs-Only-Sites-Freq-10pct.tsv"),
          path("final_output_files/${sampleId}-knownEditing-labeled.tsv"),
          path("final_output_files/${sampleId}-novelEditing.tsv"),
          emit: final_outputs

    tuple val(sampleId),
          path("summary_stats_files/${sampleId}-editingSummary.txt"),
          path("summary_stats_files/${sampleId}-FinalDiscoverySummary.txt"),
          emit: summary_stats

    // Intermediate files for publishing
    path "outTable_${sampleId}"
    path "${sampleId}-first"

    script:
    """
    echo "[Reditools Postprocessing] DNA-RNA editing output directory: ${reditoolsOutputDir}"
    echo "[Reditools Postprocessing] Starting post-processing for raw RNA editing calls for sample ${sampleId}."

    if bash reditools-v1-post-run-filtering.sh \
            ${bamFile} \
            ${reditoolsOutputDir} \
            ${genomeFa} \
            ${spliceSitesAnnot} \
            ${excludedContigs} \
            ${rmskGtf} \
            ${snpGtf} \
            ${rediportalsDbGtf} \
            ${task.cpus} \
            ${sampleId}; then
        echo "[Reditools Postprocessing] Post-processing completed successfully for sample ${sampleId}."
    else
        echo "[Reditools Postprocessing] Error: Post-processing failed for sample ${sampleId}." >&2
        exit 1
    fi
    """
}