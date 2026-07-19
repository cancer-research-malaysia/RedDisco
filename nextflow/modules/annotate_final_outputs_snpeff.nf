// Reditools v1 - VCF conversion, SnpEff annotation, and final TSV extraction
process ANNOTATE_FINAL_OUTPUTS_SNPEFF {
    tag "$sampleId"
    label 'finalAnnotateSnpeff'
    container "${params.container__snpeff}"
    publishDir { "${params.outputDir}/reditools-v1/${sampleId}/annotated-sifted-sites" }, mode: 'copy'

    input:
    tuple val(sampleId),
          path(allEditingTsv)

    output:
    // Generic intermediate VCFs mapped from the allEditing
    path "${sampleId}-allSingleSubsEditing.pre-annot.vcf"
    path "${sampleId}-allSingleSubsEditing.annot.vcf"

    // Transcript-level Specialized Subsets 
    tuple val(sampleId),
          path("${sampleId}_single-subs-all-sites.tsv"),
          path("${sampleId}_AG-subs-all-sites.tsv"),
          path("${sampleId}_AG-subs-rediportal-known-sites.tsv"),
          emit: out_tables

    script:
    def rawVcf = "${sampleId}-allSingleSubsEditing.pre-annot.vcf"
    def annVcf = "${sampleId}-allSingleSubsEditing.annot.vcf"
    """
    echo "[Step 1/3] Converting REDItools TSV to VCF..."
    bash reditools-v1-filt-output-to-custom-vcf.sh \
        ${allEditingTsv} \
        ${rawVcf}

    echo "[Step 2/3] Running SnpEff annotation..."
    # Write to a temporary file first so a crash never leaves a truncated VCF file tracked as 'valid'
    snpEff ann \
        -Xmx${task.memory.toGiga()}g \
        -v hg38 \
        ${rawVcf} > ${annVcf}.tmp
    
    mv ${annVcf}.tmp ${annVcf}

    echo "[Step 3/3] Extracting SnpEff-annotated VCF to TSV for sample ${sampleId}..."
    bash snpeff-annot-to-tsv.sh \
        ${annVcf} \
        ${sampleId}

    echo "[ANNOTATE_FINAL_OUTPUTS_SNPEFF] All steps completed successfully for sample ${sampleId}."
    """
}
