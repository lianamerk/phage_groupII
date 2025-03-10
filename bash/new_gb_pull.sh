#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <accessions_file>"
    exit 1
fi

# File containing the list of accessions
accessions_file="$1"

# Check if the file exists
if [ ! -f "$accessions_file" ]; then
    echo "File not found: $accessions_file"
    exit 1
fi

# Loop through each accession in the file and download the corresponding GenBank file
while IFS= read -r acc; do
    efetch -db nucleotide -format gb -id "$acc" > "${acc}.gb"
done < "$accessions_file"