#!/bin/bash
# Run in pharokka conda env

file_path="genomes_wth_groupII.txt"

# Loop through the file
while IFS= read -r line; do
    echo "Submitting pharokka ${line}..."
    sbatch -t 6-00:00 -p eddy -J ${line}_pharokka -c 8 -o ../data/slurms/pharokka-${line}.out --mem-per-cpu=4G --wrap="pharokka.py --force -i ../data/genomes/groupII_millard/${line}/${line}.fna -o ../data/pharokka_output/${line} -d /n/eddy_lab/data/pharokka -p ${line} -t 8"
    sleep 1
done < "$file_path"