#!/bin/bash
#SBATCH --job-name=rnaseq_pipeline
#SBATCH --output=pipeline_master_%j.log
#SBATCH --error=pipeline_master_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --partition=normal
#SBATCH --account=your_account  # Replace with your actual account

# Print job info
echo "Starting pipeline master process"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_JOB_NODELIST"
echo "Date: $(date)"

# Load modules if needed
# module load nextflow

# Or activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rnaseq-pipeline

# Set up working directory - change to your pipeline directory
cd /path/to/rnaseq-minimap-oarfish

# Create output directory for logs
mkdir -p logs

# Launch the Nextflow pipeline with cluster profile
./run.sh \
  --sample_sheet samples.csv \
  --reference reference/transcriptome.fa \
  --output results \
  --threads 16 \
  --mapq 10 \
  --params "--filter-group no-filters --model-coverage" \
  --profile cluster

# Print completion information
echo "Master process completed at: $(date)"