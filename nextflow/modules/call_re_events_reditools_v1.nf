// Defines the Reditools main process
process CALL_RE_EVENTS_REDITOOLS_V1 {
    tag "$sampleId"
    label 'callREEvents'
    container "${params.container__reditools}"
    publishDir "${params.outputDir}/reditools-v1/${sampleId}", mode: 'copy'

    input:
    // Receives the existing tuple structure
    tuple val(sampleId), path(bamFile), path(indexFile), path(dnaBamFile), path(dnaIndexFile)
    path genomeFa

    output:
    // Output the DNA-RNA editing output dir from Reditools
    tuple val(sampleId), path("DnaRna_${sampleId}"), emit: reditools_output

    script:
    """
    # run the first round of Reditools for raw RNA editing site calling
    if REDItoolDnaRna.py -i ${bamFile} -j ${dnaBamFile} -f ${genomeFa} -o ${sampleId} -t ${task.cpus} -c 1,1 -m 30,255 -v 1 -q 30,30 -e -n 0.0 -N 0.0 -u -l -p -s2 -g2 -S; then
       echo "[Reditools] First round of Reditools completed successfully for sample ${sampleId}." && DNARNA=\$(find ./${sampleId} -type d -name "DnaRna_*" | head -n1) && echo "[Reditools] DNA-RNA editing output directory: \${DNARNA}"
       mv \${DNARNA} ./DnaRna_${sampleId}

    else
       echo "[Reditools] Error: First round of Reditools failed for sample ${sampleId}." >&2
       exit 1
    fi
    """
}
