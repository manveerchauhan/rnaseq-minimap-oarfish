# RNA-seq Pipeline with Minimap2 and Oarfish

This Nextflow pipeline processes long-read RNA-seq data using Minimap2 for alignment and Oarfish for transcript quantification.

## Features

- Parallel processing of multiple RNA-seq samples
- Alignment with Minimap2 (long-read high-quality mode)
- Quality filtering (MAPQ ≥ 10) using Samtools
- SAM to BAM conversion and sorting with Samtools
- Transcript quantification with Oarfish's coverage model
- MultiQC reporting for quality assessment
- Configurable for different computing environments (local, cluster, cloud)

## Requirements

- [Nextflow](https://www.nextflow.io/) (v21.04.0 or later)
- [Minimap2](https://github.com/lh3/minimap2)
- [Samtools](http://www.htslib.org/)
- [Oarfish](https://github.com/COMBINE-lab/oarfish) (Transcript quantification tool)
- [MultiQC](https://multiqc.info/) (optional, for reports)

## Installation

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed

### Step-by-Step Installation

1. **Install Conda** (if not already installed):
   ```bash
   # Download the Miniconda installer
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
   
   # Make the installer executable
   chmod +x miniconda.sh
   
   # Run the installer
   ./miniconda.sh
   
   # Follow the prompts to complete installation
   # Then initialize conda in your shell
   source ~/.bashrc
   ```

2. **Clone this repository**:
   ```bash
   git clone https://github.com/manveerchauhan/rnaseq-minimap-oarfish.git
   cd rnaseq-minimap-oarfish
   ```

3. **Make the run script executable**:
   ```bash
   chmod +x run.sh
   ```

4. **Create and activate the conda environment**:
   ```bash
   # Create environment from the provided yml file
   conda env create -f environment.yml
   
   # Activate the environment
   conda activate rnaseq-pipeline
   ```

   The environment.yml file includes all necessary dependencies:
   - nextflow, minimap2, samtools, oarfish, multiqc (core pipeline tools)
   - Python libraries for data analysis
   - Quality control tools
   - Utilities for sequence manipulation and visualization
   
5. **Verify installation**:
   ```bash
   # Check that key tools are available
   nextflow -version
   minimap2 --version
   samtools --version
   oarfish --version
   ```

## Directory Structure

```
rnaseq-minimap-oarfish/
├── main.nf                  # Main pipeline script
├── nextflow.config          # Configuration file
├── run.sh                   # Runner script
├── environment.yml          # Conda environment file
├── data/                    # Put your FASTQ files here
│   ├── sample1_1.fastq.gz
│   └── sample1_2.fastq.gz
├── reference/               # Reference transcriptome
│   └── genome.fa
└── results/                 # Output directory
```

## Usage

### Basic Usage

```bash
./run.sh
```

### Custom Parameters

```bash
./run.sh --reads "data/*_{1,2}.fastq.gz" \
         --reference "reference/transcriptome.fa" \
         --output results \
         --threads 16 \
         --mapq 10 \
         --params "--filter-group no-filters --model-coverage" \
         --profile cluster
```

### Command-line Options

- `-h, --help`: Show help message
- `-r, --reads PATH`: Path to input reads (glob pattern)
- `-g, --reference PATH`: Path to reference transcriptome
- `-o, --output PATH`: Output directory
- `-t, --threads NUMBER`: Number of CPU threads
- `-q, --mapq NUMBER`: Minimum mapping quality score (default: 10)
- `-p, --params STRING`: Additional parameters for Oarfish
- `-c, --config PATH`: Custom Nextflow config file
- `-x, --profile STRING`: Nextflow profile (standard, cluster, cloud)

## Output

The pipeline generates the following output structure:

```
results/
├── minimap2/               # Minimap2 alignment results
│   └── sample1.sam
├── bam/                    # Filtered and sorted BAM files
│   ├── sample1.filtered.sorted.bam
│   └── sample1.filtered.sorted.bam.bai
├── oarfish/                # Oarfish transcript quantification results
│   ├── sample1.quant
│   ├── sample1.meta_info.json
│   └── sample1.ambig_info.tsv
├── multiqc/                # MultiQC reports
│   └── multiqc_report.html
└── reports/                # Nextflow execution reports
    ├── execution_report.html
    ├── timeline.html
    └── trace.txt
```

## Quick Start Guide

Once the installation is complete, follow these steps to run the pipeline:

1. **Prepare your data**:
   ```bash
   # Create directories for your data
   mkdir -p data reference results
   
   # Copy or link your FASTQ files to the data directory
   cp /path/to/your/reads/*.fastq.gz data/
   
   # Copy or link your reference transcriptome to the reference directory
   cp /path/to/your/transcriptome.fa reference/
   ```

2. **Run the pipeline**:
   ```bash
   # Basic run with default parameters
   ./run.sh
   
   # Or run with custom parameters
   ./run.sh --reads "data/*.fastq.gz" \
            --reference "reference/transcriptome.fa" \
            --output "results" \
            --threads 16 \
            --mapq 10 \
            --params "--filter-group no-filters --model-coverage"
   ```

3. **Examine the results**:
   ```bash
   # View the transcript quantification results
   head results/oarfish/sample1.quant
   
   # Open the MultiQC report in a browser
   firefox results/multiqc/multiqc_report.html
   ```

## Customization

### Using Different Profiles

The pipeline comes with predefined profiles for different environments:

- `standard`: For local execution
- `cluster`: For SLURM cluster environments
- `cloud`: For AWS Batch execution

```bash
./run.sh --profile cluster
```

### Using Containers

You can enable container support by uncommenting and modifying the container section in `nextflow.config`:

```groovy
process {
    container = 'yourusername/rnaseq-tools:latest'
}

singularity {
    enabled = true
    autoMounts = true
}
```

## Troubleshooting

### Common Issues

1. **Memory errors during execution**:
   - Adjust memory settings in `nextflow.config`
   - Try using the `-resume` flag to continue from the last successful step

2. **Missing dependencies**:
   - Ensure you've activated the conda environment: `conda activate rnaseq-pipeline`
   - Check if all tools are installed correctly: `which minimap2 samtools oarfish`

3. **Input format issues**:
   - Ensure your FASTQ files are properly formatted and not corrupted
   - Check that your reference transcriptome is in proper FASTA format

4. **Pipeline execution errors**:
   - Check the error logs in `.nextflow.log`
   - Look for specific process errors in the `work/` directory

### Getting Help

If you encounter issues not covered here, please open an issue on the GitHub repository with:
- The exact command you ran
- The complete error message
- Your system information
- Any relevant logs

## Citation

If you use this pipeline in your work, please cite:

- Minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100.
- Oarfish: Zare Jousheghani Z, Patro R (2024). Oarfish: Enhanced probabilistic modeling leads to improved accuracy in long read transcriptome quantification. bioRxiv 2024.02.28.582591; doi: https://doi.org/10.1101/2024.02.28.582591
- Nextflow: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.