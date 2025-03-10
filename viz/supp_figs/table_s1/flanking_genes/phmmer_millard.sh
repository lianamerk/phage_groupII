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
    -J "phmmer_millard_${prefix}" \
    --mem-per-cpu=4G \
    -o /dev/null \
    --wrap="phmmer -o ${prefix}_millard.out --tblout ${prefix}_millard.tblout ${faa_file} /n/eddy_lab/data/phage-2023_11/inphared_14Dec2023/14Dec2023_genomes_translated.faa"
done

