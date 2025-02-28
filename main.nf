#!/usr/bin/env nextflow

// Define parameters
params.reads = "$baseDir/data/*_{1,2}.fastq.gz"  // Input fastq files pattern
params.reference = "$baseDir/reference/genome.fa" // Reference genome
params.output_dir = "$baseDir/results"           // Output directory
params.oarfish_params = ""                        // Additional oarfish parameters
params.threads = 8                               // Number of threads
params.mapq = 10                                 // Minimum mapping quality score

log.info """\
    RNA-SEQ PIPELINE WITH MINIMAP2 AND OARFISH
    =========================================
    reads        : ${params.reads}
    reference    : ${params.reference}
    output_dir   : ${params.output_dir}
    threads      : ${params.threads}
    min_mapq     : ${params.mapq}
    """
    .stripIndent()

// Create channels for input files
Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .set { read_pairs_ch }

// Define process for alignment with minimap2 using long-read high-quality mode
process MINIMAP2_ALIGN {
    tag "$sample_id"
    publishDir "${params.output_dir}/minimap2", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads) from read_pairs_ch
    path reference from params.reference
    
    output:
    tuple val(sample_id), path("${sample_id}.sam") into alignment_ch
    
    script:
    // Handle both single and paired-end reads
    def reads_input = reads instanceof List ? "${reads[0]} ${reads[1]}" : "${reads}"
    
    """
    minimap2 -ax lr:hq --eqx -N 100 -t ${params.threads} ${reference} ${reads_input} > ${sample_id}.sam
    """
}

// Convert SAM to BAM, filter by quality, and sort
process SAM_TO_FILTERED_BAM {
    tag "$sample_id"
    publishDir "${params.output_dir}/bam", mode: 'copy'
    
    input:
    tuple val(sample_id), path(sam_file) from alignment_ch
    
    output:
    tuple val(sample_id), path("${sample_id}.filtered.sorted.bam"), path("${sample_id}.filtered.sorted.bam.bai") into sorted_bam_ch
    
    script:
    """
    # Convert SAM to BAM
    samtools view -bS ${sam_file} > ${sample_id}.bam
    
    # Filter BAM by mapping quality
    samtools view -b -q ${params.mapq} ${sample_id}.bam > ${sample_id}.filtered.bam
    
    # Sort and index the filtered BAM
    samtools sort -@ ${params.threads} ${sample_id}.filtered.bam -o ${sample_id}.filtered.sorted.bam
    samtools index ${sample_id}.filtered.sorted.bam
    """
}

// Run oarfish for RNA-seq analysis with coverage model
process OARFISH_ANALYSIS {
    tag "$sample_id"
    publishDir "${params.output_dir}/oarfish", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam_file), path(bam_index) from sorted_bam_ch
    
    output:
    path "${sample_id}*" into oarfish_results_ch
    
    script:
    """
    oarfish -j ${params.threads} -a ${bam_file} -o ${sample_id} --filter-group no-filters --model-coverage
    """
}

// Generate MultiQC report for all samples
process MULTIQC {
    publishDir "${params.output_dir}/multiqc", mode: 'copy'
    
    input:
    path '*' from oarfish_results_ch.collect()
    
    output:
    path 'multiqc_report.html'
    
    script:
    """
    multiqc ${params.output_dir}
    """
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'SUCCESS' : 'FAILED' }"
    log.info "Execution duration: $workflow.duration"
}