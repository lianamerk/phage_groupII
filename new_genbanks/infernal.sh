#!/usr/bin/env bash
#
# infernal.sh
#
# Description:
#   This script submits two Infernal cmsearch jobs to the SLURM scheduler.
#   It searches for group I and group II introns (using two different
#   covariance models) in a FASTA file of phage genomes.
#
# Usage:
#   bash infernal.sh
#
# Requirements:
#   - SLURM job scheduler (sbatch command).
#   - Infernal (cmsearch).
#   - Access to the specified covariance model (*.cm) files.
#   - Access to the specified FASTA file of genomes.

########################################
# User-adjustable variables
########################################

# SLURM parameters
JOB_TIME="6-00:00"        # 6 days
PARTITION="eddy"
JOB_CORES=64
MEM_PER_CPU="6G"
NODES=1

# Directory paths
OUT_DIR="/n/eddy_lab/users/lmerk/phage_groupII/data/infernal"
MODEL_DIR="/n/eddy_lab/users/lmerk/phage_groupII/data/infernal/models"
FASTA="/n/eddy_lab/data/phage-2023_11/inphared_14Dec2023/14Dec2023_genomes.fa"

########################################
# Job 1: Group I introns
########################################

sbatch \
  -t "${JOB_TIME}" \
  -p "${PARTITION}" \
  -J "g1_intron_millard" \
  -c "${JOB_CORES}" \
  --mem-per-cpu="${MEM_PER_CPU}" \
  -N "${NODES}" \
  -o "${OUT_DIR}/g1_intron_millard.out" \
  --wrap="cmsearch \
           -A ${OUT_DIR}/g1_intron_millard.aout \
           --tblout ${OUT_DIR}/g1_intron_millard.tblout \
           --verbose \
           -E 0.01 \
           --cpu ${JOB_CORES} \
           ${MODEL_DIR}/14x_GISSD_and_RF00028_gpI_seeds.cm \
           ${FASTA}"

########################################
# Job 2: Group II introns
########################################

sbatch \
  -t "${JOB_TIME}" \
  -p "${PARTITION}" \
  -J "g2_intron_millard" \
  -c "${JOB_CORES}" \
  --mem-per-cpu="${MEM_PER_CPU}" \
  -N "${NODES}" \
  -o "${OUT_DIR}/g2_intron_millard.out" \
  --wrap="cmsearch \
           -A ${OUT_DIR}/g2_intron_millard.aout \
           --tblout ${OUT_DIR}/g2_intron_millard.tblout \
           --verbose \
           -E 0.01 \
           --cpu ${JOB_CORES} \
           ${MODEL_DIR}/g2_introns.cm \
           ${FASTA}"
