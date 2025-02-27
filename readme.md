# RNA-seq Pipeline with Minimap2 and Oarfish

This Nextflow pipeline processes bulk RNA-seq data using Minimap2 for alignment and Oarfish for RNA-seq analysis.

## Features

- Parallel processing of multiple RNA-seq samples
- Alignment with Minimap2 (splice-aware)
- SAM to BAM conversion and sorting with Samtools
- RNA-seq analysis with Oarfish
- MultiQC reporting for quality assessment
- Configurable for different computing environments (local, cluster, cloud)

## Requirements

- [Nextflow](https://www.nextflow.io/) (v21.04.0 or later)
- [Minimap2](https://github.com/lh3/minimap2)
- [Samtools](http://www.htslib.org/)
- [Oarfish](https://github.com/oarfish) (RNA-seq analysis tool)
- [MultiQC](https://multiqc.info/) (optional, for reports)

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/manveerchauhan/rnaseq-minimap-oarfish.git
   cd rnaseq-minimap-oarfish
   ```

2. Make the run script executable:
   ```bash
   chmod +x run.sh
   ```

3. Install required tools:
   ```bash
   # Using conda (recommended)
   conda env create -f environment.yml
   conda activate rnaseq-pipeline
   
   # Or using individual installations
   conda install -c bioconda minimap2 samtools multiqc
   # Install oarfish as per their instructions
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
├── reference/               # Reference genome
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
         --reference "reference/genome.fa" \
         --output results \
         --threads 16 \
         --params "--additional-oarfish-params" \
         --profile cluster
```

### Command-line Options

- `-h, --help`: Show help message
- `-r, --reads PATH`: Path to input reads (glob pattern)
- `-g, --reference PATH`: Path to reference genome
- `-o, --output PATH`: Output directory
- `-t, --threads NUMBER`: Number of CPU threads
- `-p, --params STRING`: Additional parameters for Oarfish
- `-c, --config PATH`: Custom Nextflow config file
- `-x, --profile STRING`: Nextflow profile (standard, cluster, cloud)

## Output

The pipeline generates the following output structure:

```
results/
├── minimap2/               # Minimap2 alignment results
│   └── sample1.sam
├── bam/                    # Sorted BAM files
│   ├── sample1.sorted.bam
│   └── sample1.sorted.bam.bai
├── oarfish/                # Oarfish RNA-seq analysis results
│   └── sample1_*
├── multiqc/                # MultiQC reports
│   └── multiqc_report.html
└── reports/                # Nextflow execution reports
    ├── execution_report.html
    ├── timeline.html
    └── trace.txt
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

## Citation

If you use this pipeline in your work, please cite:

- Minimap2: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100.
- Oarfish: [Citation for Oarfish]
- Nextflow: Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316-319.
