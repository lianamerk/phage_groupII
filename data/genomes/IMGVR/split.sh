#!/bin/bash

# Define the input file and output directory
input_file="../../../IMGVR_hit_genomes.fna"
output_dir="IMGVR_hit_genomes"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Split the FASTA file
awk -v output_dir="$output_dir" '
BEGIN {
  FS="|"
}
/^>/ {
  if (seq) {
    print header "\n" seq > output_file
  }
  header = $1
  output_file = output_dir "/" substr(header, 2) ".fna"
  seq = ""
  next
}
{
  seq = seq $0 "\n"
}
END {
  if (seq) {
    print header "\n" seq > output_file
  }
}' "$input_file"

echo "Splitting completed. Individual files are saved in the '$output_dir' directory."
