#!/bin/bash
# Run in defensefinder conda env

file_path="genomes_wth_groupII.txt"

# Loop through the file
while IFS= read -r line; do
    echo "/n/eddy_lab/users/lmerk/19/bakta_output/${line}/${line}.faa"
    sbatch -t 6-00:00 -p eddy -J ${line}_b_defense -c 8 --mem-per-cpu=4G --wrap="defense-finder run ../data/bakta_output/${line}/${line}.faa -o ../data/defensefinder/bakta_defense/${line}/"
    # on bakta
    sbatch -t 6-00:00 -p eddy -J ${line}_p_defense -c 8 --mem-per-cpu=4G --wrap="defense-finder run ../data/pharokka_output/${line}/phanotate.faa -o ../data/defensefinder/pharokka_defense/${line}/"
    sleep 1
done < "$file_path"