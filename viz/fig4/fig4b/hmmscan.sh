#!/bin/bash

# Input file containing sequences
INPUT_FILE="IEP.faa"
OUTPUT_DIR="split_sequences"

# Create output directory if it does not exist
mkdir -p "$OUTPUT_DIR"

# Variables for the job submission
PARTITION="eddy"
TIME_LIMIT="6-00:00"
JOB_NAME="hmmscan"
CPU_COUNT=32
MEM_PER_CPU="6G"
PFAM_DB="/n/eddy_lab/data/pfam-35.0/Pfam-A.hmm"

# Split the multi-sequence FASTA file into individual FASTA files
awk '/^>/ {if (seq) close(out); seq++; out=sprintf("%s/seq_%d.fasta", "'"$OUTPUT_DIR"'", seq)} {print > out}' "$INPUT_FILE"

# Loop through each generated FASTA file and submit the job
for seq_file in "$OUTPUT_DIR"/*.fasta; do
    sbatch -t "$TIME_LIMIT" -p "$PARTITION" -J "$JOB_NAME" -c "$CPU_COUNT" \
           --mem-per-cpu="$MEM_PER_CPU" --wrap="hmmscan --cpu $CPU_COUNT --max --tblout ${seq_file}_pfam $PFAM_DB $seq_file"
done

echo "Job submission completed."
