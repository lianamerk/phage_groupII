#!/usr/bin/env python3
import glob
import csv
import os

# Create or overwrite the CSV file with column headers:
with open("all_headers.csv", "w", newline="") as out_csv:
    writer = csv.writer(out_csv)
    writer.writerow(["sequence_id", "file_name_no_ext"])
    
    # Loop over all FASTA files ending with "subset20.faa"
    for faa_file in glob.glob("*subset20.faa"):
        # 1) Remove the ".faa" extension 
        base_name = os.path.splitext(faa_file)[0]
        # e.g. "RVT_retron_hmmscan_subset20.faa" -> "RVT_retron_hmmscan_subset20"

        # 2) Remove the prefix "RVT_"
        if base_name.startswith("RVT_"):
            base_name = base_name.replace("RVT_", "", 1)
            # e.g. "RVT_retron_hmmscan_subset20" -> "retron_hmmscan_subset20"

        # 3) Remove "_hmmscan"
        base_name = base_name.replace("_hmmscan", "")
        # e.g. "retron_hmmscan_subset20" -> "retron_subset20"

        # 4) Remove "_subset20"
        base_name = base_name.replace("_subset20", "")
        # e.g. "retron_subset20" -> "retron"

        file_name_no_ext = base_name

        # Now read the .faa file line by line
        with open(faa_file, "r") as infile:
            for line in infile:
                if line.startswith(">"):
                    header = line[1:].strip()        # remove ">" and trailing whitespace
                    sequence_id = header.split()[0] # first token only
                    writer.writerow([sequence_id, file_name_no_ext])

print("Done! See all_headers.csv.")
