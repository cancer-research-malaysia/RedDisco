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
    path "${sampleId}-allEditing.vcf"
    path "${sampleId}-allEditing.ann.vcf"

    // Transcript-level Specialized Subsets 
    tuple val(sampleId),
          path("${sampleId}_transcript_level_editing.tsv"),
          path("${sampleId}_hyper_edited_protein-coding.tsv"),
          path("${sampleId}_isolated_high_penetrance_sites_protein-coding.tsv"),
          path("${sampleId}_neoantigen_candidates_protein-coding.tsv"),
          path("${sampleId}.known_sites_all_freq.tsv")
          emit: transcript_level_tables

    script:
    def rawVcf = "${sampleId}-allEditing.vcf"
    def annVcf = "${sampleId}-allEditing.ann.vcf"
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
