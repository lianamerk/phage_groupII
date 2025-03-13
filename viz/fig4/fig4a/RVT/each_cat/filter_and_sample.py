import os
import random
import subprocess
from Bio import SeqIO

def reformat_filter_subsample(directory):
    # Make sure "esl-reformat" is accessible (on your PATH or specify full path)
    for filename in os.listdir(directory):
        if filename.endswith(".aout"):
            # Build the input and output paths
            input_path = os.path.join(directory, filename)
            base_name = os.path.splitext(filename)[0]  # remove ".aout"
            faa_path = os.path.join(directory, f"{base_name}.faa")
            
            # 1) Convert .aout to .faa
            with open(faa_path, 'w') as faa_out:
                subprocess.run(["esl-reformat", "fasta", input_path], stdout=faa_out, check=True)
            
            # 2) Parse the new .faa
            records = list(SeqIO.parse(faa_path, "fasta"))
            
            # 3) Filter by length > 100
            filtered_records = [rec for rec in records if len(rec.seq) > 100]
            
            # Write filtered records
            filtered_outfile = os.path.join(directory, f"{base_name}_filtered.faa")
            with open(filtered_outfile, 'w') as f:
                SeqIO.write(filtered_records, f, "fasta")
            
            # 4) Take a random subsample of 20 sequences (or fewer if not enough)
            sample_size = 20 if len(filtered_records) >= 20 else len(filtered_records)
            subsample = random.sample(filtered_records, sample_size)
            
            # Write subsample
            sample_outfile = os.path.join(directory, f"{base_name}_sample20.faa")
            with open(sample_outfile, 'w') as s:
                SeqIO.write(subsample, s, "fasta")

if __name__ == "__main__":
    # Update this to your target directory containing the .aout files
    my_directory = "./"
    reformat_filter_subsample(my_directory)
