#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl = 2

// Define parameters
params.sample_sheet = "$baseDir/samples.csv"     // Sample sheet with sample IDs and FASTQ paths
params.reads = "$baseDir/data/*_{1,2}.fastq.gz"  // Legacy input fastq files pattern
params.reference = "$baseDir/reference/genome.fa" // Reference transcriptome
params.output_dir = "$baseDir/results"           // Output directory
params.oarfish_params = ""                        // Additional oarfish parameters
params.threads = 8                               // Number of threads
params.mapq = 10                                 // Minimum mapping quality score

log.info """\
    RNA-SEQ PIPELINE WITH MINIMAP2 AND OARFISH
    =========================================
    sample_sheet : ${params.sample_sheet}
    reference    : ${params.reference}
    output_dir   : ${params.output_dir}
    threads      : ${params.threads}
    min_mapq     : ${params.mapq}
    """
    .stripIndent()

// Define process for alignment with minimap2 using long-read high-quality mode and streaming directly to BAM
process MINIMAP2_ALIGN {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(reads)
    path reference
    
    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai")
    
    script:
    // Handle both single and paired-end reads
    def reads_input = reads instanceof List ? "${reads[0]} ${reads[1]}" : "${reads}"
    
    """
    # Align and directly convert to sorted BAM in one pipeline
    minimap2 -ax lr:hq --eqx -N 100 -t ${params.threads} ${reference} ${reads_input} | \
    samtools view -bS - | \
    samtools sort -@ ${params.threads} -o ${sample_id}.sorted.bam -
    
    # Index the BAM file
    samtools index ${sample_id}.sorted.bam
    """
}

// Filter BAM by quality and prepare for analysis
process SAM_TO_FILTERED_BAM {
    tag "$sample_id"
    publishDir "${params.output_dir}/bam", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index)
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered.sorted.bam"), path("${sample_id}.filtered.sorted.bam.bai")
    
    script:
    """
    # Filter BAM by mapping quality
    samtools view -b -q ${params.mapq} ${bam_file} > ${sample_id}.filtered.sorted.bam
    
    # Index the filtered BAM
    samtools index ${sample_id}.filtered.sorted.bam
    """
}

// Run oarfish for RNA-seq analysis with coverage model
process OARFISH_ANALYSIS {
    tag "$sample_id"
    publishDir "${params.output_dir}/oarfish", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index)
    
    output:
    path "${sample_id}*"
    
    script:
    """
    oarfish -j ${params.threads} -a ${bam_file} -o ${sample_id} --filter-group no-filters --model-coverage ${params.oarfish_params}
    """
}

// Define the workflow
workflow {
    // Create input channel from sample sheet if it exists
    if (file(params.sample_sheet).exists()) {
        samples_ch = Channel
            .fromPath(params.sample_sheet)
            .splitCsv(header: true)
            .map { row -> 
                def sample_id = row.sample_id
                def fastq_files = row.fastq_path.split(',')
                
                if (fastq_files.size() == 1) {
                    // Single-end reads
                    return [sample_id, file(fastq_files[0])]
                } else if (fastq_files.size() >= 2) {
                    // Paired-end reads - take first two files
                    return [sample_id, [file(fastq_files[0]), file(fastq_files[1])]]
                } else {
                    error "Invalid entry in sample sheet for sample: ${sample_id}"
                }
            }
        
        log.info "Using sample sheet: ${params.sample_sheet}"
    } else {
        // Fallback to legacy mode with glob pattern
        log.info "Sample sheet not found, using legacy input pattern: ${params.reads}"
        
        samples_ch = Channel
            .fromFilePairs(params.reads, checkIfExists: true)
    }

    // Reference as channel
    reference_ch = Channel.value(file(params.reference))
    
    // Run the pipeline
    alignment_ch = MINIMAP2_ALIGN(samples_ch, reference_ch)
    bam_ch = SAM_TO_FILTERED_BAM(alignment_ch)
    OARFISH_ANALYSIS(bam_ch)
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
    log.info "Execution duration: $workflow.duration"
}