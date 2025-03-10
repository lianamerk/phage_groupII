#!/usr/bin/env bash

# Loop over all .faa files in the current directory
for faa_file in *.faa
do
  # Remove the .faa extension to get the base filename
  prefix="${faa_file%.faa}"
  
  # Submit the job
  sbatch \
    -p eddy \
    -t 60 \
    -c 12 \
    -n 1 \
    -J "phmmer_${prefix}" \
    --mem-per-cpu=4G \
    -o /dev/null \
    --wrap="phmmer -o ${prefix}.out --tblout ${prefix}.tblout ${faa_file} /n/eddy_lab/data/uniprot-2024nov/uniprot.fa"
done

