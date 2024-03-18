# Destination directory where the copied genomes will be saved
fasta_directory="/n/eddy_lab/users/lmerk/genomes/groupII_millard/"

# File containing the list
file_path="genomes_wth_groupII.txt"

# Loop through the file
while IFS= read -r line; do
    echo "${fasta_directory}${line}/${line}.fna"
    sbatch -t 6-00:00 -p eddy -J ${line} -c 64 --mem-per-cpu=6G --wrap="esl-sfetch --index ${fasta_directory}${line}/${line}.fna"
    sleep 1
done < "$file_path"

