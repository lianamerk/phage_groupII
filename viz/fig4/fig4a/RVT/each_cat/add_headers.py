#!/usr/bin/env python3

import glob
import os

def prepend_filename_to_headers():
    for faa_file in glob.glob("*hmmscan_sample20.faa"):
        # e.g. "RVT_g2l_hmmscan_sample20.faa" -> "RVT_g2l_hmmscan_sample20"
        base_name = os.path.splitext(faa_file)[0]
        
        # Remove the known prefixes and suffixes
        # "RVT_g2l_hmmscan_sample20" -> "g2l"
        prefix = base_name.replace("RVT_", "", 1).replace("_hmmscan_sample20", "")
        
        # Output file name, e.g. "RVT_g2l_hmmscan_sample20_prefixed.faa"
        out_file = f"{base_name}_prefixed.faa"
        
        with open(faa_file, "r") as infile, open(out_file, "w") as outfile:
            for line in infile:
                if line.startswith(">"):
                    # Prepend "g2l_" (example) to the existing header (minus the '>')
                    outfile.write(">" + prefix + "_" + line[1:])
                else:
                    outfile.write(line)
        
        print(f"Processed {faa_file} -> {out_file}")

if __name__ == "__main__":
    prepend_filename_to_headers()
