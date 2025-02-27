#!/bin/bash

# Run the RNA-seq Nextflow pipeline

# Default values
READS="data/*_{1,2}.fastq.gz"
REFERENCE="reference/genome.fa"
OUTPUT_DIR="results"
THREADS=8
OARFISH_PARAMS=""

# Display help message
function show_help {
  echo "Usage: ./run.sh [options]"
  echo "RNA-seq pipeline with Minimap2 and Oarfish"
  echo ""
  echo "Options:"
  echo "  -h, --help                Show this help message"
  echo "  -r, --reads PATH          Path to input reads (default: data/*_{1,2}.fastq.gz)"
  echo "  -g, --reference PATH      Path to reference genome (default: reference/genome.fa)"
  echo "  -o, --output PATH         Output directory (default: results)"
  echo "  -t, --threads NUMBER      Number of threads (default: 8)"
  echo "  -p, --params STRING       Additional oarfish parameters"
  echo "  -c, --config PATH         Custom config file"
  echo "  -x, --profile STRING      Nextflow profile (standard, cluster, cloud)"
  echo ""
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -h|--help)
      show_help
      exit 0
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
      OUTPUT_DIR="$2"
      shift 2
      ;;
    -t|--threads)
      THREADS="$2"
      shift 2
      ;;
    -p|--params)
      OARFISH_PARAMS="$2"
      shift 2
      ;;
    -c|--config)
      CONFIG="-c $2"
      shift 2
      ;;
    -x|--profile)
      PROFILE="-profile $2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      show_help
      exit 1
      ;;
  esac
done

# Create directories if they don't exist
mkdir -p ${OUTPUT_DIR}

# Run the Nextflow pipeline
echo "Starting RNA-seq pipeline with Minimap2 and Oarfish..."
nextflow run main.nf \
  --reads "${READS}" \
  --reference "${REFERENCE}" \
  --output_dir "${OUTPUT_DIR}" \
  --threads ${THREADS} \
  --oarfish_params "${OARFISH_PARAMS}" \
  ${CONFIG} ${PROFILE} \
  -resume

echo "Pipeline execution completed."
