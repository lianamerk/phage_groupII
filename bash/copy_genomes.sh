# Destination directory where the copied genomes will be saved
destination_directory="/n/eddy_lab/users/lmerk/genomes/groupII_millard/"

list_file="./genomes_wth_groupII.txt"


# Loop through the list file and copy each genome directory
while IFS= read -r genome_dir || [[ -n "$genome_dir" ]]; do
    cp -r "/n/eddy_lab/data/IMGVR-2022-12-19_7$genome_dir" "$destination_directory"
done < "$list_file"