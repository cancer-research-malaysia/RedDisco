// Defines the Reditools main process
process CALL_RE_EVENTS_REDITOOLS_V1 {
    tag "$sample_id"
    label 'call-re-events-reditools'
    container "${params.container__reditools}"
    publishDir "${params.outputDir}/reditools-v1/${sample_id}", mode: 'copy'

    input:
    // Receives the existing tuple structure
    tuple val(sample_id), path(bam_file), path(bai_file)

    output:
    // Output the HTML/ZIP files
    tuple val(sample_id), path("*fastqc.zip"), path("*fastqc.html")

    script:
    """
    # run first round of reditools
    """
}