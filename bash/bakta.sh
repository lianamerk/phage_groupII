#!/bin/bash
# Run in bakta conda env

# File containing the list
file_path="genomes_wth_groupII.txt"

# Loop through the file
while IFS= read -r line; do
    echo "Submitting bakta ${line}..."
    sbatch -t 6-00:00 -p eddy -J ${line} -c 64 --mem-per-cpu=6G -o ../data/slurms/bakta-${line}.out --wrap="bakta --force --db /n/eddy_lab/users/lmerk/databases/full/db -o ../data/bakta_output/${line} ../data/genomes/groupII_millard/${line}/${line}.fna"
    sleep 1
done < "$file_path"
