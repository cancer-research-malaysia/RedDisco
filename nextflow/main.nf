#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import original modules
include { CALL_RE_EVENTS_REDITOOLS_V1 } from './modules/call_re_events_reditools_v1.nf'
include { POSTPROCESS_REDITOOLS_V1_OUTPUTS } from './modules/postprocess_reditools_v1_outputs.nf'
include { ANNOTATE_FINAL_OUTPUTS_SNPEFF } from './modules/annotate_final_outputs_snpeff.nf'

// Import lifted modules (with FastQC aliased for raw and trimmed stages)
include { TRIM_READS_FASTP } from './modules/trim_reads_fastp.nf'
include { QUANT_COUNTS_KALLISTO } from './modules/quant_counts_kallisto.nf'
include { QC_READS_FASTQC as QC_RAW } from './modules/qc_reads_fastqc.nf'
include { QC_READS_FASTQC as QC_TRIMMED } from './modules/qc_reads_fastqc.nf'


workflow {

    // Input handling for BAM files from metadata
    if (!params.manifestPath) {
        error "--manifestPath must be provided"
    }
    // Validate genome reference exists
    if (!params.genomeFa) {
        error "Genome reference file (--genomeFa) must be provided"
    }

    // Channel 1: Extract BAM information
    bam_ch = channel
        .fromPath(params.manifestPath)
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.bamPath }
        .map { row ->
            tuple(
                row.sampleName,
                file(row.bamPath),
                file(row.baiPath),
                file(row.dnaBamPath),
                file(row.dnaBaiPath)
            )
        }
    
    // Channel 2: Extract FASTQ information + Strict Naming/Pairing Checks
    fastq_ch = channel
        .fromPath(params.manifestPath)
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.fastq1 && row.fastq2 } 
        .map { row ->
            def read1 = file(row.fastq1)
            def read2 = file(row.fastq2)
            
            // Re-apply the strict R1/R2 suffix and extension validation from the directory loop
            def r1 = read1.name =~ /[_\.][Rr]?1\.(fastq|fq)(\.gz)?$/
            def r2 = read2.name =~ /[_\.][Rr]?2\.(fastq|fq)(\.gz)?$/
            
            if (!r1 || !r2) {
                error "Invalid read pairing format in manifest for sample ${row.sampleName}: ${read1.name}, ${read2.name}"
            }
            
            tuple(row.sampleName, read1, read2)
        }

    // Reference files setup
    genome_fa                = file(params.genomeFa)
    splice_sites             = file(params.spliceSitesAnnotation)
    excluded_contigs         = file(params.excludedContigs)
    rmsk_gtf                 = file(params.rmskGtf)
    rmsk_gtf_index           = file(params.rmskGtfIndex)
    snp_gtf                  = file(params.snpGtf)
    snp_gtf_index            = file(params.snpGtfIndex)
    rediportals_db_gtf       = file(params.rediportalsDbGtf)
    rediportals_db_gtf_index = file(params.rediportalsDbGtfIndex)

    ///////////////////////// --- Main Analysis --- /////////////////////////////

    // --- TRACK A: Trimming, QC, & Kallisto Pseudoalignment ---
    
    // 1. Run FastQC on Raw Reads
    QC_RAW(fastq_ch, "raw")

    // 2. Run Trimming
    TRIM_READS_FASTP(fastq_ch)

    // 3. Run FastQC on Trimmed Reads
    QC_TRIMMED(TRIM_READS_FASTP.out.trimmed_reads, "trimmed")

    // 4. Run Kallisto quantification (Always runs if data is present)
    if (!params.kallistoIndex) { 
        error "Missing mandatory parameter: --kallistoIndex required for transcriptome pseudoalignment." 
    }
    kallisto_index_ch = channel.value(file(params.kallistoIndex))
    
    QUANT_COUNTS_KALLISTO(
        TRIM_READS_FASTP.out.trimmed_reads,
        kallisto_index_ch
    )


    // --- TRACK B: Main REDItools RNA-Editing Variant Pipeline ---

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
            .map { sampleId, allEditing, _agSubs, _knownLabeled, _novel ->
                   tuple(sampleId, allEditing) }
    )

    // Completion handler
    workflow.onComplete = {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
        log.info "Duration: $workflow.duration"
    }
}
