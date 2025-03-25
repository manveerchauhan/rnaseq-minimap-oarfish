#!/bin/bash

# Default parameter values
READS=""
SAMPLE_SHEET=""
REFERENCE="$PWD/reference/genome.fa"
OUTPUT="$PWD/results"
THREADS=8
MAPQ=10
PARAMS=""
CONFIG=""
PROFILE="standard"

# Display help
function show_help {
    echo "RNA-Seq Pipeline with Minimap2 and Oarfish"
    echo ""
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -h, --help                 Show this help message"
    echo "  -s, --sample_sheet PATH    Path to sample sheet CSV file"
    echo "  -r, --reads PATH           Path to input reads (glob pattern) - legacy mode"
    echo "  -g, --reference PATH       Path to reference transcriptome"
    echo "  -o, --output PATH          Output directory"
    echo "  -t, --threads NUMBER       Number of CPU threads"
    echo "  -q, --mapq NUMBER          Minimum mapping quality score"
    echo "  -p, --params STRING        Additional parameters for Oarfish"
    echo "  -c, --config PATH          Custom Nextflow config file"
    echo "  -x, --profile STRING       Nextflow profile (standard, cluster, cloud)"
    echo ""
    echo "Example using sample sheet:"
    echo "  $0 --sample_sheet samples.csv --reference reference/transcriptome.fa --threads 16"
    echo ""
    echo "Example using legacy mode:"
    echo "  $0 --reads \"data/*.fastq.gz\" --reference reference/transcriptome.fa --threads 16"
    exit 1
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -h|--help)
            show_help
            ;;
        -s|--sample_sheet)
            SAMPLE_SHEET="$2"
            shift 2
            ;;
        -r|--reads)
            READS="$2"
            shift 2
            ;;
        -g|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -q|--mapq)
            MAPQ="$2"
            shift 2
            ;;
        -p|--params)
            PARAMS="$2"
            shift 2
            ;;
        -c|--config)
            CONFIG="$2"
            shift 2
            ;;
        -x|--profile)
            PROFILE="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            ;;
    esac
done

# Validate inputs
if [[ -z "$SAMPLE_SHEET" && -z "$READS" ]]; then
    echo "ERROR: Either --sample_sheet or --reads must be provided"
    show_help
fi

# Build Nextflow command
NF_CMD="nextflow run main.nf"

# Add parameters
if [[ ! -z "$SAMPLE_SHEET" ]]; then
    NF_CMD="$NF_CMD --sample_sheet $SAMPLE_SHEET"
fi

if [[ ! -z "$READS" ]]; then
    NF_CMD="$NF_CMD --reads '$READS'"
fi

NF_CMD="$NF_CMD --reference $REFERENCE"
NF_CMD="$NF_CMD --output_dir $OUTPUT"
NF_CMD="$NF_CMD --threads $THREADS"
NF_CMD="$NF_CMD --mapq $MAPQ"

if [[ ! -z "$PARAMS" ]]; then
    NF_CMD="$NF_CMD --oarfish_params '$PARAMS'"
fi

if [[ ! -z "$CONFIG" ]]; then
    NF_CMD="$NF_CMD -c $CONFIG"
fi

NF_CMD="$NF_CMD -profile $PROFILE"

# Print the command
echo "Executing: $NF_CMD"

# Execute the command
eval $NF_CMD