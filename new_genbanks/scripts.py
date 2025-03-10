import numpy as np
import pandas as pd
from Bio import SeqIO, GenBank
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import subprocess

def extract_hit_df_easy(filename, filt = True):
    initial = 0
    colspecs = []

    with open(filename, 'r') as infile:
        data = infile.readlines()
        for i in data[1:2]:
            widths = [len(col)+1 for col in i.split()]
            # print(i.split())

    for n in range(0, len(widths)):
        if n == len(widths)-1:
            colspecs.append((initial, -1))
        else:
            colspecs.append((initial, initial + widths[n]))
        initial +=  widths[n]

    hit_df = pd.read_fwf(filename, skipfooter=10, skiprows=2, colspecs = colspecs, names=['target_name', 'accession', 'query_name', 'accession_2', 'mdl', \
                                                        'mdl_from', 'mdl_to', 'seq_from', 'seq_to', 'strand', 'trunc', 'pass', 'gc', 'bias', \
                                                        'score', 'E-value', 'inc', 'description_of_target'])
    if filt == True:
        hit_df = hit_df.loc[hit_df['E-value'] < 0.01]
        
    # Some are on +, some -, make a stard and end to make comparing easier later
    hit_df['start'] = hit_df[['seq_from','seq_to']].min(axis=1)
    hit_df['end'] = hit_df[['seq_from','seq_to']].max(axis=1)

    hit_df['length'] = hit_df['end'] - hit_df['start']
    
    return hit_df
            
def dereplicate_hits(df):
    derepped_hits = pd.DataFrame()
    df = df.reset_index(drop=True).rename_axis('hit').reset_index()
    
    targets = list(df['target_name'].unique())

    index_to_keep = []
    index_to_ignore = []

    while len(targets) > 0:
        target = targets.pop(0)
        sub_df = df.loc[df['target_name'] == target]                        
        # For each hit
        for index, row in sub_df.iterrows():
            # Set the current hit to be the current row
            curr_row = row
            # For each hit
            for test_index, test_row in sub_df.iterrows():

                # That isn't the current hit we are testing against
                if index == test_index:
                    pass
                elif test_index in index_to_ignore or index in index_to_ignore:
                    pass
                else:
                    # if the one we are testing start is within the region of the hit we chose in the outer loop
                    if (row['start']-200 < test_row['start'] < row['end']+200) or (row['start'] - 200 < test_row['end'] < row['end']+200):  

                        # and if the test e value is worse (ie greater than the hit we chose in the outer loop)
                        if test_row['E-value'] > row['E-value']:
                            # Drop the test hit from the dataframe
                            # sub_df = sub_df.drop(test_index)
                            index_to_keep.append(index)
                            index_to_ignore.append(test_index)

                            
                        else:
                            # If the row e-value is higher, drop that one instead
                            # sub_df = sub_df.drop(index, errors='ignore')
                            index_to_keep.append(test_index)
                            index_to_ignore.append(index)
                            
            index_to_ignore = list(set(index_to_ignore))

    df = df.sort_values('target_name')
    return df.loc[~df['hit'].isin(index_to_ignore)], df.loc[df['hit'].isin(index_to_ignore)]


def collapse_hits(g2_df):
    rows_to_add = []
    index_to_drop = []

    for target in g2_df.target_name.unique():
        sub_df = g2_df.loc[g2_df['target_name'] == target]   

        window_size = 3000

        for index, row in sub_df.iterrows():
            current_value = row['seq_from']
            indices_within_window = sub_df[(sub_df['seq_from'] >= current_value - window_size) & (sub_df['seq_from'] <= current_value + window_size)].index

            if len(indices_within_window) == 2:
                row_to_add = g2_df.iloc[indices_within_window[0]].copy()

                low_bound = np.min([g2_df.iloc[indices_within_window[0]].seq_from, g2_df.iloc[indices_within_window[1]].seq_from])
                high_bound = np.max([g2_df.iloc[indices_within_window[0]].seq_to, g2_df.iloc[indices_within_window[1]].seq_to])

                row_to_add['seq_from'] = low_bound
                row_to_add['seq_to'] = high_bound
                row_to_add['query_name'] = 'collapsed_g2'

                # Remove the original rows
                index_to_drop.append(indices_within_window[0])
                index_to_drop.append(indices_within_window[1])

                rows_to_add.append(row_to_add)

    g2_df = g2_df.drop(list(set(index_to_drop)))
    g2_df = g2_df.append(rows_to_add, ignore_index=True).drop_duplicates()
    g2_df = g2_df.sort_values('target_name')
    
    return g2_df


def update_genbank_only_g2(genbank_directory, all_genomes, intron_df):
    """
    File will be saved as "genbank_file_modified"
    """

    # Iterate over grouped hits
    for genome in all_genomes:
        # Load GenBank file for the current genome
        record = SeqIO.read(genome, 'genbank')
        
        genome_name = genome.split('/')[1].split('.gb')[0]

        curr_introns = intron_df.loc[intron_df.target_name == genome_name]
        
        # Iterate over intron hits in the group
        for index, row in curr_introns.iterrows():
            # Ensure start and end positions are in the correct order
            start = int(row['seq_from'])
            end = int(row['seq_to'])
            if start > end:
                start, end = end, start

            # Create a single feature spanning from start to end position
            feature = SeqFeature(FeatureLocation(start, end), type='ncRNA')

            # Add qualifiers to the feature
            feature.qualifiers['label'] = [row['query_name']]
            feature.qualifiers['query_name'] = [row['accession_2']]
            feature.qualifiers['E-value'] = [row['E-value']]

            if row['strand'] =='-':
                feature.strand = -1
            else:
                feature.strand = +1
                
            # Append the new feature to the existing features list in the GenBank record
            record.features.append(feature)

        # Save the modified GenBank record
        modified_genbank_file = os.path.join(f'genbank_with_intron_hits/{genome_name}_with_intron.gb')
        SeqIO.write(record, modified_genbank_file, 'genbank')
        