#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { CALL_RE_EVENTS_REDITOOLS_V1 } from './modules/call_re_events_reditools_v1.nf'
include { POSTPROCESS_REDITOOLS_V1_OUTPUTS } from './modules/postprocess_reditools_v1_outputs.nf'
include { ANNOTATE_FINAL_OUTPUTS_SNPEFF } from './modules/annotate_final_outputs_snpeff.nf'

workflow {

    // Input handling for BAM files from metadata
    if (!params.manifestPath) {
        error "--manifestPath must be provided"
    }

    // Read from manifest file with BAM information
    bam_ch = channel
        .fromPath(params.manifestPath)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.sampleName,
                file(row.bamPath),
                file(row.baiPath),
                file(row.dnaBamPath),
                file(row.dnaBaiPath)
            )
        }

    // Validate genome reference exists
    if (!params.genomeFa) {
        error "Genome reference file (--genomeFa) must be provided"
    }

    genome_fa                = file(params.genomeFa)
    splice_sites             = file(params.spliceSitesAnnotation)
    excluded_contigs         = file(params.excludedContigs)
    rmsk_gtf                 = file(params.rmskGtf)
    rmsk_gtf_index           = file(params.rmskGtfIndex)
    snp_gtf                  = file(params.snpGtf)
    snp_gtf_index            = file(params.snpGtfIndex)
    rediportals_db_gtf       = file(params.rediportalsDbGtf)
    rediportals_db_gtf_index = file(params.rediportalsDbGtfIndex)

    // --- Main Analysis ---

    // Step 1: Call RNA editing events using REDItools v1
    CALL_RE_EVENTS_REDITOOLS_V1(
        bam_ch,
        genome_fa
    )

    // Step 2: Join bam_ch and reditools output on sampleId before postprocessing
    // This guarantees correct sample pairing regardless of parallel execution order
    postprocess_input_ch = bam_ch
        .join(CALL_RE_EVENTS_REDITOOLS_V1.out.reditools_output, by: 0)

    POSTPROCESS_REDITOOLS_V1_OUTPUTS(
        postprocess_input_ch.map { sampleId, bam, bai, dnaBam, dnaBai, _reditoolsOut ->
            tuple(sampleId, bam, bai, dnaBam, dnaBai) },
        postprocess_input_ch.map { sampleId, _bam, _bai, _dnaBam, _dnaBai, reditoolsOut ->
            tuple(sampleId, reditoolsOut) },
        genome_fa,
        splice_sites,
        excluded_contigs,
        rmsk_gtf,
        rmsk_gtf_index,
        snp_gtf,
        snp_gtf_index,
        rediportals_db_gtf,
        rediportals_db_gtf_index
    )

    // Step 3: VCF conversion, SnpEff annotation, and final TSV extraction
    // Single-emit source — no join needed
    ANNOTATE_FINAL_OUTPUTS_SNPEFF(
        POSTPROCESS_REDITOOLS_V1_OUTPUTS.out.final_outputs
            .map { sampleId, _allEditing, agSubs, _knownLabeled, _novel ->
                   tuple(sampleId, agSubs) }
    )

    // Completion handler
    workflow.onComplete = {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
        log.info "Duration: $workflow.duration"
    }
}
