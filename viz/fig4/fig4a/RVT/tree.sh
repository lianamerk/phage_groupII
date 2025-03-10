#!/bin/bash

#SBATCH -t 6-00:00:00          # Job will run for up to 1 day
#SBATCH -p eddy                # Partition name
#SBATCH -J mafft_iqtree        # Job name
#SBATCH -c 64                  # Number of cores
#SBATCH --mem-per-cpu=6G       # Memory per CPU core
#SBATCH --output=mafft_iqtree_%j.out  # Output file for job (stdout)
#SBATCH --error=mafft_iqtree_%j.err   # Error file for job (stderr)

# Step 1: Align the sequences using MAFFT
mafft --thread 64 --auto each_cat/both.faa > aligned_both.faa

# # Step 2: Build the phylogenetic tree using IQ-TREE
iqtree2 -s aligned_both.faa -m MFP -bb 1000 -nt 64 -pre aligned_both
