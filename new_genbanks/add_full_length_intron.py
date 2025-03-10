import os
import glob
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import pandas as pd

# Step 1: Read the CSV file
csv_path = 'g2_df.csv'  # Update with the actual path
data = pd.read_csv(csv_path)
output_dir = "genbank_with_full_intron"
os.makedirs(output_dir, exist_ok=True)
data = data.dropna(subset=["Accession"])

def add_group2_feature(genbank_file, accession, subtype, start, stop):
    records = list(SeqIO.parse(genbank_file, "genbank"))

    # Strand determination
    strand_var = 1
    if start > stop:
        start, stop = stop, start
        strand_var = -1
    
    # Create the feature
    group2_feature = SeqFeature(
        FeatureLocation(start, stop),
        type="group2",
        strand=strand_var,
        qualifiers={
            "subtype": subtype,
        },
    )
    
    # Append to the first record
    records[0].features.append(group2_feature)
    
    # Output file
    output_file = os.path.join(output_dir, f"{accession}_full_intron.gb")
    SeqIO.write(records, output_file, "genbank")
    print(f"Saved: {output_file}")

for index, row in data.iterrows():
    accession = row["Accession"]
    pattern = f'genbank_with_intron_hits/{accession}*.gb*'
    matching_files = glob.glob(pattern)
    
    if not matching_files:
        print(f"No GenBank file found for accession {accession} (pattern: {pattern})")
        continue
    
    # Usually, there's only one matching file, but if there's more, you could loop or pick the first
    file = matching_files[0]
    
    # Safely parse start/end
    try:
        start = int(row["Start"])
        stop = int(row["End"])
    except ValueError:
        print(f"Invalid Start/End for {accession}. Skipping.")
        continue
    
    # Add the feature
    try:
        add_group2_feature(file, accession, row["Subtype"], start, stop)
    except Exception as e:
        print(f"Error adding group2 feature for {accession}: {e}")
        continue