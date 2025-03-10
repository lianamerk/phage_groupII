#!/usr/bin/env python3

# File: filter_fasta_by_prefix.py
import sys

if len(sys.argv) < 4:
    print("Usage: python filter_fasta_by_prefix.py <fasta_in> <prefixes_file> <fasta_out>")
    sys.exit(1)

fasta_in = sys.argv[1]
prefixes_file = sys.argv[2]
fasta_out = sys.argv[3]

# Read all the prefixes into a set for quick lookup
prefixes = set()
with open(prefixes_file) as pf:
    for line in pf:
        line = line.strip()
        if line:
            prefixes.add(line)

keep = False
with open(fasta_in) as fin, open(fasta_out, "w") as fout:
    for line in fin:
        # If this is a header line (starts with ">")
        if line.startswith(">"):
            # Extract the header ID (e.g. "fig|1179224.4.peg.4192")
            # Usually the ID is everything after ">" up to the first whitespace
            header = line[1:].rstrip()  # no trailing newline
            if any(pref in header for pref in prefixes):
                # keep the record

                keep = True
                fout.write(line)  # write the header
            else:
                keep = False
        else:
            # If we're keeping this record, write the sequence lines as well
            if keep:
                fout.write(line)
