# RNA-seq Pipeline with Minimap2 and Oarfish

A Nextflow pipeline for long-read RNA-seq data processing using Minimap2 for alignment and Oarfish for transcript quantification.

## Features

- Nextflow DSL2 workflow for improved modularity
- Processes multiple RNA-seq samples in parallel
- Alignment with Minimap2 (long-read high-quality mode)
- Quality filtering and BAM conversion with Samtools
- Transcript quantification with Oarfish
- Compatible with both compressed (`.fastq.gz`) and uncompressed (`.fastq`) files

## Quick Start

```bash
# Clone the repository
git clone <repository-url>
cd rnaseq-minimap-oarfish

# Make run script executable
chmod +x run.sh

# Activate conda environment
conda activate bulk-rnaseq-longbench

# Run the pipeline
./run.sh --sample_sheet samples.csv \
         --reference /path/to/reference/transcriptome.fa \
         --output results \
         --threads 16 \
         --mapq 10 \
         --params "--filter-group no-filters --model-coverage" \
         --profile standard
```

## Sample Sheet Format

Create a CSV file with your samples:

```
sample_id,fastq_path
sample1,/path/to/sample1.fastq.gz
sample2,/path/to/sample2.fastq
sample3,/path/to/sample3_1.fastq.gz,/path/to/sample3_2.fastq.gz
```

## Running on HPC

1. Create a SLURM submission script:

```bash
#!/bin/bash
#SBATCH --partition="partition_name"
#SBATCH --nodes=1
#SBATCH --account="your_account"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

# Load modules
module load Anaconda3

# Activate conda environment
conda activate bulk-rnaseq-longbench

# Run the pipeline
bash /path/to/run.sh --sample_sheet samples.csv \
     --reference /path/to/reference.fa \
     --output results \
     --threads 16 \
     --mapq 10 \
     --params "--filter-group no-filters --model-coverage" \
     --profile standard
```

2. Submit the job:

```bash
sbatch submit_job.slurm
```

## Resource Requirements

- **CPU**: 8-16 cores recommended
- **Memory**: 32-64GB for medium datasets, 100GB+ for large datasets
- **Disk**: 5-10x the size of your input data
- **Time**: Varies by sample size; 6-24 hours per sample for alignment and quantification

## Pipeline Outputs

```
results/
├── minimap2/              # Minimap2 alignment results (.sam)
├── bam/                   # Filtered and sorted BAM files
│   ├── sample1.filtered.sorted.bam
│   └── sample1.filtered.sorted.bam.bai
└── oarfish/               # Oarfish transcript quantification results
    ├── sample1.quant
    ├── sample1.meta_info.json
    └── sample1.ambig_info.tsv
```

## Command-line Options

- `--sample_sheet`: Path to the sample sheet CSV
- `--reference`: Path to reference transcriptome
- `--output`: Output directory
- `--threads`: Number of CPU threads
- `--mapq`: Minimum mapping quality score
- `--params`: Additional parameters for Oarfish
- `--profile`: Execution profile (standard, cluster)

## Troubleshooting

- **Permission denied**: Use `chmod +x run.sh` or run with `bash run.sh`
- **Nextflow version error**: This pipeline uses DSL2; ensure Nextflow 22.04.0 or later
- **Memory issues**: Increase memory allocation in your job submission
- **Missing files**: Check your sample sheet paths; use absolute paths for reliability

## Dependencies

- Nextflow (v22.04.0 or later)
- Minimap2
- Samtools
- Oarfish

These can be installed via the included conda environment file.