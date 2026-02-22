// Reditools v1 - VCF conversion, SnpEff annotation, and final TSV extraction
process ANNOTATE_FINAL_OUTPUTS_SNPEFF {
    tag "$sampleId"
    label 'finalAnnotateSnpeff'
    container "${params.container__snpeff}"
    publishDir "${params.outputDir}/reditools-v1/${sampleId}/final-annotated-SnpEff", mode: 'copy'

    input:
    tuple val(sampleId),
          path(agSubsOnlyTsv)

    output:
    // Intermediate VCF files
    path "${sampleId}-AG-Subs-Only-Sites-Freq-10pct.vcf"
    path "${sampleId}-AG-Subs-Only-Sites-Freq-10pct.ann.vcf"

    // Site-level TSV outputs
    tuple val(sampleId),
          path("${sampleId}.known_sites_freq10pct.tsv"),
          path("${sampleId}.known_sites_freq10pct_genic.tsv"),
          emit: site_level_tables

    // Gene-level TSV outputs
    tuple val(sampleId),
          path("${sampleId}_gene_level_editing.tsv"),
          path("${sampleId}_gene_level_editing_subset_protein_coding_gt2_sites.tsv"),
          emit: gene_level_tables

    script:
    def rawVcf = "${sampleId}-AG-Subs-Only-Sites-Freq-10pct.vcf"
    def annVcf = "${sampleId}-AG-Subs-Only-Sites-Freq-10pct.ann.vcf"
    """
    echo "[Step 1/3] Converting REDItools TSV to VCF with strand correction for sample ${sampleId}..."
    bash reditools-v1-filt-output-to-custom-vcf.sh \
        ${agSubsOnlyTsv} \
        ${rawVcf}

    echo "[Step 2/3] Running SnpEff annotation for sample ${sampleId}..."
    snpEff ann \
        -Xmx${task.memory.toGiga()}g \
        -v hg38 \
        ${rawVcf} > ${annVcf}

    echo "[Step 3/3] Extracting SnpEff-annotated VCF to TSV for sample ${sampleId}..."
    bash snpeff-annot-to-tsv.sh \
        ${annVcf} \
        ${sampleId}

    echo "[ANNOTATE_FINAL_OUTPUTS_SNPEFF] All steps completed successfully for sample ${sampleId}."
    """
}
