#!/bin/bash

# SLURM parameters
PARTITION="eddy"
TIME_LIMIT="6-00:00" # 6 days
CPUS=64
MEM_PER_CPU="6G"
TARGET_FILE="/n/eddy_lab/data/phage-2023_11/inphared_14Dec2023/14Dec2023_genomes.fa"

# Directories
FASTA_DIR="./"
OUTPUT_DIR="./"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Iterate over each .fasta file in ../
for fasta_file in "$FASTA_DIR"*.fasta; do
  # Get the base name of the file (no path or extension)
  fasta_basename=$(basename "$fasta_file")

  # Output file naming
  output_file="$OUTPUT_DIR/${fasta_basename}.nhmmer.out"
  log_file="$OUTPUT_DIR/${fasta_basename}.out"

  # Construct SLURM sbatch command
  sbatch \
    -t "$TIME_LIMIT" \
    -p "$PARTITION" \
    -J "nhmmer_${fasta_basename}" \
    -c "$CPUS" \
    --mem-per-cpu="$MEM_PER_CPU" \
    -N 1 \
    -o "$log_file" \
    --wrap="nhmmer --tblout $output_file $fasta_file $TARGET_FILE"

  echo "Submitted job for $fasta_file"
done

echo "All jobs submitted."

