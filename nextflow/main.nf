#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { CALL_RE_EVENTS_REDITOOLS_V1 } from './modules/call_re_events_reditools_v1.nf'
include { POSTPROCESS_REDITOOLS_V1_OUTPUTS } from './modules/postprocess_reditools_v1_outputs.nf'
include { WRANGLE_FINAL_OUTPUTS_REDITOOLS_V1 } from './modules/wrangle_final_outputs_reditools_v1.nf'

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
    genome_fa = file(params.genomeFa)
    splice_sites = file(params.spliceSitesAnnotation)
    excluded_contigs = file(params.excludedContigs)
    rmsk_gtf = file(params.rmskGtf)
    rmsk_gtf_index = file(params.rmskGtfIndex)
    snp_gtf = file(params.snpGtf)
    snp_gtf_index = file(params.snpGtfIndex)
    rediportals_db_gtf = file(params.rediportalsDbGtf)
    rediportals_db_gtf_index = file(params.rediportalsDbGtfIndex)

    // --- Main Analysis ---
    
    // Call RNA editing events using Reditools
    CALL_RE_EVENTS_REDITOOLS_V1(
        bam_ch,
        genome_fa
    )

    // Post-process Reditools outputs
    POSTPROCESS_REDITOOLS_V1_OUTPUTS(
        bam_ch,
        CALL_RE_EVENTS_REDITOOLS_V1.out.reditools_output,
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
    // Combine the outputs using join on the common sampleId (the first element)
    wrangle_input = POSTPROCESS_REDITOOLS_V1_OUTPUTS.out.for_wrangling
        .join(POSTPROCESS_REDITOOLS_V1_OUTPUTS.out.out_tables)

    // Use an explicit map with named parameters instead of 'it'
    WRANGLE_FINAL_OUTPUTS_REDITOOLS_V1(
        wrangle_input.map { sampleId, known, pos, pos_alu, _firstalu, _second -> 
            tuple(sampleId, known, pos, pos_alu) 
        },
        wrangle_input.map { sampleId, _known, _pos, _pos_alu, firstalu, second -> 
            tuple(sampleId, firstalu, second) 
        }
    )

    // Completion handler
    workflow.onComplete = {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
        log.info "Duration: $workflow.duration"
    }
}
