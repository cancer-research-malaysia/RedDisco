#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { TRIM_READS_FASTP } from './modules/trim_reads_fastp.nf'
// include { QC_READS_FASTQC } from './modules/qc_reads_fastqc.nf'
// include { QUANT_COUNTS_KALLISTO } from './modules/quant_counts_kallisto.nf'

workflow {
    
    // Input handling
    if (params.manifestPath) {
        // Read from manifest file
        fastq_ch = channel
            .fromPath(params.manifestPath)
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                tuple(row.sampleName, file(row.read1Path), file(row.read2Path))
            }
    } else if (params.inputDir) {
        // Read from directory - handle multiple naming conventions
        // Supports: _R1/_R2, _r1/_r2, _1/_2 (all case-insensitive)
        // Supports: .fastq.gz, .fq.gz, .fastq, .fq
        fastq_ch = channel
            .fromFilePairs("${params.inputDir}/*_{R,r,}{1,2}.{fastq,fq}{.gz,}", size: 2, flat: true)
            .map { sampleName, read1, read2 ->
                // Verify we have proper R1/R2 pairing
                def r1 = read1.name =~ /[_\.][Rr]?1\.(fastq|fq)(\.gz)?$/
                def r2 = read2.name =~ /[_\.][Rr]?2\.(fastq|fq)(\.gz)?$/
                
                if (!r1 || !r2) {
                    error "Invalid read pairing for sample ${sampleName}: ${read1.name}, ${read2.name}"
                }
                
                tuple(sampleName, read1, read2)
            }.view()
    } else {
        error "Either --inputDir or --manifestPath must be provided"
    }

    // --- Main Analysis ---
    
    //

    // Completion handler
    workflow.onComplete = {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
        log.info "Duration: $workflow.duration"
    }
}
