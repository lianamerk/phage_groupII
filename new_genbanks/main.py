import os
import glob
import pandas as pd
import numpy as np

# Biopython
from Bio import SeqIO, GenBank
from Bio.SeqFeature import SeqFeature, FeatureLocation

# Local scripts
import scripts
import importlib

# Reload in case you've made changes in 'scripts.py'
importlib.reload(scripts)

# --- Define directories ---
data_dir = '/n/eddy_lab/users/lmerk/phage_groupII/data/'
infernal_dir = os.path.join(data_dir, 'infernal/')
genomes_dir = os.path.join(data_dir, 'genomes/')
metadata_dir = os.path.join(genomes_dir, 'inphared_metadata/')
script_output = os.path.join(data_dir, 'script_output/')
genbank_out_directory = os.path.join(script_output, "updated_genomes")

# Create output directory if it doesn't exist
os.makedirs(genbank_out_directory, exist_ok=True)

# --- Filenames ---
g1_intron_hits = os.path.join(infernal_dir, 'g1_intron_millard.tblout')
g2_intron_hits = os.path.join(infernal_dir, 'g2_intron_millard.tblout')
metadata_path = os.path.join(metadata_dir, '14Dec2023_data.tsv')
actual_genomes_path = os.path.join(genomes_dir, 'actual_genomes.txt')

# --- Load data ---
metadata = pd.read_csv(metadata_path, sep='\t').rename(columns={'Accession': 'target_name'})
actual_genomes = pd.read_csv(actual_genomes_path, header=None)[0].unique()

# Infernal file output
g1_df_unfiltered_csv = os.path.join(script_output, 'g1_df_unfiltered.csv')
g1_hits_fail_csv = os.path.join(script_output, 'g1_hits_fail.csv')
g1_hits_pass_csv = os.path.join(script_output, 'g1_hits_pass.csv')

g2_df_unfiltered_csv = os.path.join(script_output, 'g2_df_unfiltered.csv')
g2_hits_pass_csv = os.path.join(script_output, 'g2_hits_pass.csv')

# Flags / parameters
preprocess = True

# --- Preprocessing block ---
if preprocess:
    # 1. Parse group I intron hits
    g1_hits = scripts.extract_hit_df_easy(g1_intron_hits)
    g1_df_unfiltered = pd.merge(g1_hits, metadata, on='target_name', how='left')

    # 2. Parse group II intron hits
    g2_hits = scripts.extract_hit_df_easy(g2_intron_hits)
    g2_df_unfiltered = pd.merge(g2_hits, metadata, on='target_name', how='left')

    print('Dereplicating group I...')
    g1_hits_pass, g1_hits_fail = scripts.dereplicate_hits(g1_df_unfiltered)

    print('Saving group I...')
    g1_df_unfiltered.sort_values('target_name').to_csv(g1_df_unfiltered_csv, index=False)
    g1_hits_fail.to_csv(g1_hits_fail_csv, index=False)

    g1_hits_pass = g1_hits_pass.sort_values(by=['target_name', 'seq_from'])
    # Create 'intronID'
    g1_hits_pass['intronID'] = g1_hits_pass['target_name'] + '_I_' + g1_hits_pass['seq_from'].astype(str)
    # Reorder columns
    column_order = ['intronID'] + [col for col in g1_hits_pass.columns if col != 'intronID']
    g1_hits_pass = g1_hits_pass[column_order]
    # Only keep introns from actual_genomes
    g1_hits_pass = g1_hits_pass.loc[g1_hits_pass['target_name'].isin(actual_genomes)]
    g1_hits_pass.to_csv(g1_hits_pass_csv, index=False)

    print('Collapsing group II...')
    g2_df = g2_df_unfiltered.copy()
    g2_hits_pass = scripts.collapse_hits(g2_df)

    print('Saving group II...')
    g2_hits_pass = g2_hits_pass.sort_values(by=['target_name', 'seq_from'])
    g2_hits_pass['intronID'] = g2_hits_pass['target_name'] + '_II_' + g2_hits_pass['seq_from'].astype(str)
    column_order = ['intronID'] + [col for col in g2_hits_pass.columns if col != 'intronID']
    g2_hits_pass = g2_hits_pass[column_order]
    g2_hits_pass = g2_hits_pass.loc[g2_hits_pass['target_name'].isin(actual_genomes)]
    g2_hits_pass.to_csv(g2_hits_pass_csv, index=False)

    # Save all group II hits (unfiltered)
    g2_df_unfiltered = g2_df_unfiltered.sort_values(by=['target_name', 'seq_from'])
    g2_df_unfiltered['intronID'] = g2_df_unfiltered['target_name'] + '_II_' + g2_df_unfiltered['seq_from'].astype(str)
    column_order = ['intronID'] + [col for col in g2_df_unfiltered.columns if col != 'intronID']
    g2_df_unfiltered = g2_df_unfiltered[column_order]
    g2_df_unfiltered.to_csv(g2_df_unfiltered_csv, index=False)

    print('Done!')

else:
    # If not preprocessing, just read the existing CSVs
    g1_df_unfiltered = pd.read_csv(g1_df_unfiltered_csv)
    g2_df_unfiltered = pd.read_csv(g2_df_unfiltered_csv)
    g1_hits_pass = pd.read_csv(g1_hits_pass_csv)
    g2_hits_pass = pd.read_csv(g2_hits_pass_csv)


# --- Generate updated GenBank files ---
all_genomes = glob.glob('genbank_0/*.gb')
scripts.update_genbank_only_g2('genbank_0', all_genomes, g2_df_unfiltered)       