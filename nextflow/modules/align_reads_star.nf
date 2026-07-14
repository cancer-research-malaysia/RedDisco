//
process ALIGN_READS_STAR {
    // maxForks 10
    cpus params.numCores
    
    label 'alignReadsStar'
    container "${params.container__preproc}"
    publishDir { "${params.outputDir}/STAR-alignments/${sampleName}/processed_calls" }, mode: 'copy'
    
    input:
        tuple val(sampleName), path(trimmedReads)
        path starIndex

    output:
        tuple val(sampleName), path("${sampleName}-STAR_Aligned.coordSorted.bam"), path("${sampleName}-STAR_Aligned.coordSorted.bam.bai"), emit: final_bam

    script:
    """
    # variables
    SAMPLE_ID=${sampleName}
    CORES=${task.cpus}
    STAR_INDEX=${starIndex}
    READ1=${trimmedReads[0]}  # First file in the nested list will be read 1 file
    READ2=${trimmedReads[1]}

    echo "Processing files of sample \${SAMPLE_ID}"
    echo "Number of cores to use: ${task.cpus}"
    echo "The index path: \${STAR_INDEX}"
    
    # STAR normal alignment for general tools
    if star-align-REDi.sh "\${READ1}" "\${READ2}" "\${SAMPLE_ID}" ${task.cpus} "\${STAR_INDEX}"; then
        echo "STAR general alignment is complete!"
        # index only as the STAR command outputs bam sorted by coordinates already - create files with the names Nextflow expects
        samtools index -@ ${task.cpus} "\${SAMPLE_ID}-STAR*.bam"
    else
        echo "STAR alignment failed. Check logs. Exiting..."
        exit 1
    fi

    """
}
