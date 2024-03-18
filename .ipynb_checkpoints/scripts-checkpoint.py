import numpy as np
import pandas as pd
from Bio import SeqIO, GenBank
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import holoviews as hv
import iqplot
import subprocess


def extract_hit_df(filename, filt = True):
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
    
    hit_df['desc'] = hit_df['description_of_target'].str.split(',').str[0]
    hit_df['desc'] = hit_df['description_of_target'].str.split(' ').str[1:].str.join(" ")
    hit_df['desc'] = hit_df['desc'].str.replace(' genome assembly', '')
    hit_df['desc'] = hit_df['desc'].str.replace('MAG: ', '')
    hit_df['desc'] = hit_df['desc'].str.replace(' genomic sequence.', '')
    hit_df['desc'] = hit_df['desc'].str.replace(' complete genome.', '')
    hit_df['desc'] = hit_df['desc'].str.replace(' phage', '')
    hit_df['desc'] = hit_df['desc'].str.replace(' chromosome: I.', '')
    hit_df['desc'] = hit_df['desc'].str.replace(' partial genome.', '')
    hit_df['desc'] = hit_df['desc'].str.replace(' sequence.', '')
    hit_df['desc'] = hit_df['desc'].str.replace('complete sequence.', '')
    hit_df['desc'] = hit_df['desc'].str.replace('Uncultured', 'CrAssphage')


    hit_df['desc'] = hit_df['desc'].str.replace(',', '')
    hit_df['host'] = hit_df['desc'].str.split(' ').str[0]
    hit_df['name'] = hit_df['desc'].str.split(' ').str[1]

    
    return hit_df


def get_surrounding_window(intron_df, genbank_dir):
    """
    Extract a window from a FASTA file.

    Parameters:
    - intron_df (dataframe): like g2_hits_pass, the cleaned output of extract hits
    - genbank_dir (string): directory of genome directories, that then contain {sequence}.fa

    Note you need to wait a few seconds after running (check your jobs) to make sure all the
    genome windows were extracted.
    """
        
    for target_name in intron_df['target_name'].unique():
        target_dir = os.path.join('./gII_intron_6kb_windows', target_name)
        os.makedirs(target_dir, exist_ok=True)

    sbatch_params = "-p eddy -t 60 -c 4 -n 1 --mem-per-cpu=6G"

    # Iterate through rows in the DataFrame
    for index, row in intron_df.iterrows():
        fasta_file = os.path.join(genbank_dir, f'{row.target_name}/{row.target_name}.fna')

        # need to play with window so you are only looking at encoded vs what its inserted in
        window = 6000

        # Extract window from nucleotide sequence
        intron_df.at[index, 'surrounding_window'] = extract_window(fasta_file, row.seq_from - window, row.seq_to + window, row.strand)

        # Create a directory for each target_name (if not already created)
        target_dir = os.path.join('./gII_intron_6kb_windows', row.target_name)

        os.makedirs(target_dir, exist_ok=True)

        # Write the nucleotide sequence to a temporary FASTA file
        temp_nucleotide_fasta = os.path.join(target_dir, f'{row.intronID}_intron_window.fna')

        with open(temp_nucleotide_fasta, "w") as fasta_file:
            fasta_file.write(f">{row.target_name}\n{intron_df.at[index, 'surrounding_window']}")

        # Translate nucleotide sequence to protein sequence using esl-translate
        # -l 
        temp_protein_fasta = os.path.join(target_dir, f'{row.intronID}_intron_window.faa')
        translate_command = f"esl-translate {temp_nucleotide_fasta} > {temp_protein_fasta}"
        sbatch_command = f"sbatch {sbatch_params} --wrap='{translate_command}'"
        subprocess.run(sbatch_command, shell=True)
        
def hmmscan_a_window(intron_df):
    """
    hmmscan a window from a FASTA file against Pfam-A

    Parameters:
    - intron_df (dataframe): like g2_hits_pass, the cleaned output of extract hits

    Note you need to wait a few seconds after running (check your jobs).
    """
    
    sbatch_params = "-p eddy -t 60 -c 4 -n 1 --mem-per-cpu=6G"

    
    for index, row in intron_df.iterrows():
        target_dir = os.path.join('./gII_intron_6kb_windows', row.target_name)

        temp_protein_fasta = os.path.join(target_dir, f'{row.intronID}_intron_window.faa')

        hmmscan_script = f"hmmscan -o {os.path.join(target_dir, f'{row.intronID}_hmmscan_output.out')} --domtblout {os.path.join(target_dir, f'{row.intronID}_hmmscan_output.domtblout')} --tblout {os.path.join(target_dir, f'{row.intronID}_hmmscan_output.tblout')} /n/eddy_lab/data/pfam-35.0/Pfam-A.hmm {temp_protein_fasta}"
        sbatch_command = f"sbatch {sbatch_params} --wrap='{hmmscan_script}'"
        subprocess.run(sbatch_command, shell=True)
        


def hmmscan_window(intron_df, genbank_dir, window_dir, window=5000):
    """
    This combines the above 2 functions efficiently
    Extract a window from a FASTA file and hmmscan it
    Be sure to index before (see the esl-index bash script)

    Parameters:
    - intron_df (dataframe): like g2_hits_pass, the cleaned output of extract hits
    - genbank_dir (string): directory of genome directories, that then contain {sequence}.fa

    Note you need to wait a few seconds after running (check your jobs) to make sure all the
    genome windows were extracted.
    """
        
    for target_name in intron_df['target_name'].unique():
        target_dir = os.path.join(window_dir, target_name)
        os.makedirs(target_dir, exist_ok=True)

    # don't save the slurm outputs with -o /dev/null
    # sbatch_params = "-p eddy -t 60 -c 4 -n 1 --mem-per-cpu=6G"
    sbatch_params = "-p eddy -t 60 -c 4 -n 1 --mem-per-cpu=6G -o /dev/null"


    # Iterate through rows in the DataFrame
    for index, row in intron_df.iterrows():
        fasta_file = os.path.join(genbank_dir, f'{row.target_name}/{row.target_name}.fna')

        # Create a directory for each target_name (if not already created)
        target_dir = os.path.join(window_dir, row.target_name)
        os.makedirs(target_dir, exist_ok=True)
        faa_file = os.path.join(target_dir, f'{row.intronID}.faa')

        start = row.start-window
        end = row.end+window

        if start - 1 < 1:
            start = 1
        if end > row['Genome Length (bp)']:
            end = row['Genome Length (bp)']

        cmd = f"esl-sfetch -c {start}..{end} {fasta_file} {row.target_name} | esl-translate -l 30 - | tee {faa_file} | hmmscan -o {os.path.join(target_dir, f'{row.intronID}_hmmscan_output.out')} --domtblout {os.path.join(target_dir, f'{row.intronID}_hmmscan_output.domtblout')} --tblout {os.path.join(target_dir, f'{row.intronID}_hmmscan_output.tblout')} /n/eddy_lab/data/pfam-35.0/Pfam-A.hmm -"
        sbatch_command = f"sbatch {sbatch_params} --wrap='{cmd}'"
        subprocess.run(sbatch_command, shell=True)
            
            
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
            # print(f'NOW STARTING {index}')
            # For each hit
            for test_index, test_row in sub_df.iterrows():

                # That isn't the current hit we are testing against
                if index == test_index:
                    pass
                elif test_index in index_to_ignore or index in index_to_ignore:
                    pass
                else:
                    # if index == 265 and test_index == 64:
                    # print(row['start'] < test_row['start'] < row['end'])

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
                            

                        # print(index, row['start'], row['end'], row['E-value'], 'AGAINST', \
                        #       test_index, test_row['start'], test_row['end'], test_row['E-value'], \
                        #      test_row['E-value'] < row['E-value'], index_to_keep)
            # print('SUBDF is:', len(sub_df), '\n\n\n\n')

            # print(list(set(index_to_keep)))
            # derepped_hits = derepped_hits.append(sub_df)
            index_to_ignore = list(set(index_to_ignore))
            

    df = df.sort_values('target_name')
    return df.loc[~df['hit'].isin(index_to_ignore)], df.loc[df['hit'].isin(index_to_ignore)]


def get_genes_from_genbank(hit_df_pass, genbank_dir):
    description = []
    all_names = hit_df_pass['target_name'].unique()

    for name in all_names:
        play = hit_df_pass.loc[hit_df_pass['target_name'] == name]

        cds = []
        genome_dir = os.path.join(genbank_dir, name, f'{name}.gbk')

        records = SeqIO.parse(genome_dir, "genbank")

        for record in records:
            for feature in record.features:
                if feature.type == 'CDS':
                    cds.append(feature)

            for feature in cds:
                for intron_start in play.seq_from.unique():
                    if intron_start - 1000 <= feature.location.start <= intron_start + 1000:
                        desc = feature.qualifiers['product'][0]
                        if 'reductase' in desc:
                            desc = 'reductase'
                        elif 'nuclease' in desc:
                            desc = 'nuclease'
                        elif 'tail' in desc:
                            desc = 'assembly'
                        elif 'capsid' in desc:
                            desc = 'assembly'
                        elif 'baseplate' in desc:
                            desc = 'assembly'
                        elif 'head' in desc:
                            desc = 'assembly'
                        elif 'polymerase' in desc:
                            desc = 'polymerase'
                        elif 'portal' in desc:
                            desc = 'assembly'
                        elif 'structural' in desc:
                            desc = 'assembly'
                        elif 'prophage' in desc:
                            desc = 'prophage-related'
                        elif 'hydrolase' in desc:
                            desc = 'hydrolase'
                        elif 'helicase' in desc:
                            desc = 'helicase'

                        description.append(desc)
                for intron_end in play.seq_to.unique():
                    if intron_end - 1000 <= feature.location.start <= intron_end + 1000:
                        desc = feature.qualifiers['product'][0]
                        if 'reductase' in desc:
                            desc = 'reductase'
                        elif 'nuclease' in desc:
                            desc = 'nuclease'
                        elif 'tail' in desc:
                            desc = 'assembly'
                        elif 'capsid' in desc:
                            desc = 'assembly'
                        elif 'baseplate' in desc:
                            desc = 'assembly'
                        elif 'head' in desc:
                            desc = 'assembly'
                        elif 'polymerase' in desc:
                            desc = 'polymerase'
                        elif 'portal' in desc:
                            desc = 'assembly'
                        elif 'structural' in desc:
                            desc = 'assembly'
                        elif 'prophage' in desc:
                            desc = 'prophage-related'
                        elif 'hydrolase' in desc:
                            desc = 'hydrolase'
                        elif 'helicase' in desc:
                            desc = 'helicase'
                        description.append(desc)

    return description


def plot_histogram(hit_df, intron_type):
    grouped = hit_df.groupby('target_name').count().reset_index().sort_values('query_name', ascending=False)[['query_name', 'target_name']]
    data = grouped['query_name'].values

    # p = iqplot.histogram(data, bins=len(list(set(data)))+10, rug=False, title='Number of Group I Introns Per Genome')
    p = iqplot.spike(data,
                     title=f'Number of Group {intron_type} Introns Per Genome',
                     x_axis_label = 'Introns Per Genome',
                     y_axis_label = 'Number of Genomes',
                     width = 600,
                     style="spike",
                     # y_axis_type = 'log',
                    line_kwargs=dict(line_width = 20, color = 'purple',))

    tick_vals = np.arange(np.floor(data.min()), np.ceil(data.max()) + 1)
    p.xaxis.ticker = tick_vals

    # Format the x-axis ticks as integers
    p.xaxis.major_label_overrides = {tick: str(int(tick)) for tick in tick_vals}

    return p


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


def extract_window(fasta_file, start, end, strand):
    """
    Extract a window from a FASTA file.

    Parameters:
    - fasta_file (str): Path to the FASTA file.
    - start (int): Start position of the window (1-based index).
    - end (int): End position of the window (1-based index).

    Returns:
    - str: Extracted sequence window.
    """
    with open(fasta_file, 'r') as fasta_handle:
        for record in SeqIO.parse(fasta_handle, 'fasta'):
            
            if start - 1 < 1:
                start = 1
            if end > len(record.seq):
                end = len(record.seq)
            if strand == '-':
                # Extract the reverse complement for the negative strand
                window_sequence = record.seq[start - 1:end].reverse_complement()
            else:
                window_sequence = record.seq[start - 1:end]
            return str(window_sequence)
    
    
    

        
        
def map_coordinates(protein_coords, genome_start, genome_end, strand):
    mapped_start = genome_start + protein_coords[0] - 1
    mapped_end = genome_start + protein_coords[1] - 1

    # look if the start is greater than the stop, if so mark it as a negative strand
    if strand == "-":
        mapped_start, mapped_end = genome_end - protein_coords[0] + 1, genome_end - protein_coords[1] + 1
    
    return [mapped_start, mapped_end]


def parse_fasta(filename):
    protein_coords = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('>'):
                header_parts = line.strip()[1:].split(' ')
                protein_name = header_parts[0]
                coords_string = next(part for part in header_parts if part.startswith('coords=')).split('=')[1]
                start, end = map(int, coords_string.split('..'))
                protein_coords[protein_name] = [start, end]
    return protein_coords


def update_genbank_only_g2(genbank_directory, features_df, intron_df):
    """
    File will be saved as "genbank_file_modified"
    """
    
    # Group hits by genome
    grouped_hits = features_df.groupby('genome')

    # Iterate over grouped hits
    for genome, hits_group in grouped_hits:
        # Load GenBank file for the current genome
        genbank_file = os.path.join(genbank_directory, f"{genome}/{genome}.gbff")
        record = SeqIO.read(genbank_file, 'genbank')

        # Iterate over hits in the group
        for index, row in hits_group.iterrows():
            # Ensure start and end positions are in the correct order
            start = int(row['genome_start'])
            end = int(row['genome_end'])
            if start > end:
                start, end = end, start

            # Create a single feature spanning from start to end position
            feature = SeqFeature(FeatureLocation(start, end), type='other')

            # Add qualifiers to the feature
            feature.qualifiers['label'] = [row['target_name']]
            feature.qualifiers['description_of_target'] = [row['description_of_target']]
            feature.qualifiers['E-value'] = [row['E-value']]
            if row['strand'] =='-':
                feature.strand = -1
            else:
                feature.strand = +1

            # Append the new feature to the existing features list in the GenBank record
            record.features.append(feature)

        curr_introns = intron_df.loc[intron_df.target_name == genome]
        
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
        modified_genbank_file = os.path.join(genbank_directory, genome, f"{genome}_modified.gbff")
        SeqIO.write(record, modified_genbank_file, 'genbank')
    
def update_genbank(genbank_directory, features_df, g1_intron_df, g2_intron_df, defense_df):
    """
    File will be saved as "genbank_file_modified"
    """
    
    # if 2 intron hits are next to each other, a given orf is in feature table twice. remove these
    features_df = features_df.drop_duplicates(subset=['target_name', 'genome_start', 'genome_end'])

    
    # Group hits by genome
    grouped_hits = features_df.groupby('genome')
    # Iterate over grouped hits
    for genome, hits_group in grouped_hits:
        # print(genome)
        # Load GenBank file for the current genome
        genbank_file = os.path.join(genbank_directory, f"{genome}_double_annotated.gbff")
        record = SeqIO.read(genbank_file, 'genbank')

        # Iterate over hits in the group
        for index, row in hits_group.iterrows():
            # Ensure start and end positions are in the correct order
            start = int(row['genome_start'])
            end = int(row['genome_end'])
            if start > end:
                start, end = end, start

            if row['source'] == '.tblout':
                # Create a single feature spanning from start to end position
                feature = SeqFeature(FeatureLocation(start, end), type='gene')
                feature.qualifiers['label'] = [row['query_name']]
            else:
                feature = SeqFeature(FeatureLocation(start, end), type='unsure')
                feature.qualifiers['label'] = [row['target_name']]
                feature.qualifiers['query_name'] = [row['query_name']]
                

            # Add qualifiers to the feature
            # feature.qualifiers['label'] = [row['target_name']]
            # print(row)
            # feature.qualifiers['query_name'] = [row['query_name']]
            feature.qualifiers['description_of_target'] = [row['description_of_target']]
            feature.qualifiers['E-value'] = [row['E-value']]
            feature.qualifiers['source'] = [row['source']]
            if row['strand'] =='-':
                feature.strand = -1
            else:
                feature.strand = +1

            # Append the new feature to the existing features list in the GenBank record
            record.features.append(feature)

        curr_defense = defense_df.loc[defense_df.genome == genome]

        # Iterate over intron hits in the group
        for index, row in curr_defense.iterrows():
            # Ensure start and end positions are in the correct order
            start = int(row['hit_start_nt'])
            end = int(row['hit_end_nt'])
            if start > end:
                start, end = end, start

            # Create a single feature spanning from start to end position
            feature = SeqFeature(FeatureLocation(start, end), type='gene')

            # Add qualifiers to the feature
            feature.qualifiers['label'] = [row['gene_name']]
            feature.qualifiers['i_eval'] = [row['i_eval']]

            if row['strand'] =='-':
                feature.strand = -1
            else:
                feature.strand = +1

            # Append the new feature to the existing features list in the GenBank record
            record.features.append(feature)



        curr_introns = g2_intron_df.loc[g2_intron_df.target_name == genome]

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

        curr_introns = g1_intron_df.loc[g1_intron_df.target_name == genome]

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
        modified_genbank_file = os.path.join(genbank_directory, f"{genome}_double_and_intron_and_defense.gbff")
        SeqIO.write(record, modified_genbank_file, 'genbank')
    
def hit_to_genome(intron_df):
    all_hit_df = pd.DataFrame()

    window = 6000
    for index, row in intron_df.iterrows():
        target_dir = os.path.join('./gII_intron_6kb_windows', row.target_name)

        window_start = row.start - window
        window_end = row.end + window

        strand = row.strand

        hit_path = os.path.join(target_dir, f'{row.intronID}_hmmscan_output.tblout')
        prot_fasta = os.path.join(target_dir, f'{row.intronID}_intron_window.faa')
        parsed_fasta = parse_fasta(prot_fasta)


        initial = 0
        colspecs = []

        with open(hit_path, 'r') as infile:
            data = infile.readlines()
            for i in data[2:3]:
                widths = [len(col)+1 for col in i.split()]

        for n in range(0, len(widths)):
            if n == len(widths)-1:
                colspecs.append((initial, -1))
            else:
                colspecs.append((initial, initial + widths[n]))
            initial +=  widths[n]

        hit_df = pd.read_fwf(hit_path, skipfooter=10, skiprows=3, colspecs = colspecs, names=['target_name', 'accession', 'query_name', 'accession_2', 'E-value', \
                                                            'score', 'bias', 'E-value_1_domain', 'score_1_domain', 'bias_1_domain', 'exp_domain', 'reg_domain', \
                                                              'clu_domain', 'ov_domain', 'env_domain', 'dom_domain', 'rep_domain', 'inc_domain', 'description_of_target'])

        if len(hit_df.target_name) < 0:
            print(row.target_name)
        else:
            starts = []
            ends = []
            for i, hit in hit_df.iterrows():
                curr = map_coordinates(parsed_fasta[hit.query_name], window_start, window_end, strand)
                starts.append(curr[0])
                ends.append(curr[1])

            hit_df['genome_start'] = starts
            hit_df['genome_end'] = ends

            hit_df['strand'] = [strand for i in range(len(hit_df))]
            hit_df['genome'] = [row.target_name for i in range(len(hit_df))]

            filt_hit_df = hit_df[hit_df['E-value'] < 0.01]

            # print(filt_hit_df)
            best_rows = filt_hit_df.groupby('query_name').apply(lambda x: x.nsmallest(3, 'E-value'))[['genome', 'target_name', 'query_name', 'E-value', 'description_of_target', 'genome_start', 'genome_end', 'strand']]

            all_hit_df = all_hit_df.append(best_rows, ignore_index=True)

    all_hit_df['description_of_target'] = all_hit_df['description_of_target'].str.split(n=1).str[1]
    
    return all_hit_df

def parse_fasta_fixed(filename):
    protein_coords = {}
    strands = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('>'):
                header_parts = line.strip()[1:].split(' ')
                protein_name = header_parts[0]
                coords_string = next(part for part in header_parts if part.startswith('coords=')).split('=')[1]
                start, end = map(int, coords_string.split('..'))
                
                
                genome_loc = next(part for part in header_parts if part.startswith('source')).split('/')[1]
                genome_start = int(genome_loc.split('-')[0])
                genome_end = int(genome_loc.split('-')[1])
                
                frame = int(next(part for part in header_parts if part.startswith('frame=')).split('=')[1])
                if frame < 4:
                    strands[protein_name] = '+'
                else:
                    strands[protein_name] = '-'
                    
                actual_start = genome_start + start
                actual_end = genome_start + end

                protein_coords[protein_name] = [actual_start, actual_end]
    return (protein_coords, strands)


def hit_to_genome_with_domtbl(intron_df, outdir, genomes_dir):
    all_hit_df = pd.DataFrame()
    all_dom_df = pd.DataFrame()

    for index, row in intron_df.iterrows():
        strand = row.strand
        target_dir = os.path.join(outdir, row.target_name)
        prot_fasta = os.path.join(target_dir, f'{row.intronID}.faa')
        parsed_fasta, strands = parse_fasta_fixed(prot_fasta)

        hit_path = os.path.join(target_dir, f'{row.intronID}_hmmscan_output.tblout')
        domtbl_path = os.path.join(target_dir, f'{row.intronID}_hmmscan_output.domtblout')


        # Determine the maximum width of any line in the file
        max_line_width = max(len(line.strip()) for line in open(domtbl_path))

        # Read the first 3 lines of the file to extract column widths
        with open(domtbl_path, 'r') as infile:
            for i, line in enumerate(infile):
                if i == 2:  # We're interested in the third line (header)
                    # Split the line into columns and calculate the width for each column
                    widths = [len(col) for col in line.split()]  # Widths of each column in the header
                    total_width = sum(widths) + len(widths) - 1  # Total width including spaces between columns
                    break  # Exit loop after processing the header

        # Define colspecs based on extracted header widths
        start = 0
        colspecs = []
        for width in widths:
            colspecs.append((start, start + width))
            start += width + 1  # Add 1 for the space between columns

        # Extend the width of the last column to the maximum line width
        colspecs[-1] = (colspecs[-1][0], max_line_width)

        # Read the fixed-width file into a DataFrame
        domtbl_df = pd.read_fwf(domtbl_path, skipfooter=10, skiprows=3, colspecs=colspecs, 
                                 names=['target_name', 'accession', 'tlen', 'query_name', 'accession_2', 
                                        'qlen', 'E-value', 'score_full', 'bias_full', '#', 'of', 'c-Evalue', 
                                        'i-Evalue', 'score_1_domain', 'bias_1_domain', 'from_hmm', 'to_hmm', 
                                        'from_ali', 'to_ali', 'from_env', 'to_env', 'acc', 'description_of_target'])


        if len(domtbl_df.target_name) < 0:
            print(row.target_name)
        else:
            starts = []
            ends = []
            strand_list = []
            for i, hit in domtbl_df.iterrows():
                curr = parsed_fasta[hit.query_name]
                curr_strand = strands[hit.query_name]

                strand_list.append(curr_strand)
                nt_start = hit.from_ali * 3
                nt_end = hit.to_ali * 3

                # if hit.accession == 'PF00078.30':
                #     print(curr, curr[1] - nt_start, curr[1] - nt_end)

                if curr_strand == '+':
                    starts.append(curr[0] + nt_start)
                    ends.append(curr[0]+ nt_end)
                else:
                    starts.append(curr[0] - nt_start)
                    ends.append(curr[0] - nt_end)

            domtbl_df['genome_start'] = starts
            domtbl_df['genome_end'] = ends

            domtbl_df['strand'] = strand_list
            domtbl_df['genome'] = [row.target_name for i in range(len(domtbl_df))]

            filt_dom_df = domtbl_df[domtbl_df['E-value'] < 0.01]
            best_rows = filt_dom_df.groupby('query_name').apply(lambda x: x.nsmallest(3, 'E-value'))[['genome', 'accession', 'target_name', 'query_name', 'E-value', 'description_of_target', 'genome_start', 'genome_end', 'strand']]
            best_rows['intronID'] = row.intronID
            all_dom_df = all_dom_df.append(best_rows, ignore_index=True)

        initial = 0
        colspecs = []

        with open(hit_path, 'r') as infile:
            data = infile.readlines()
            for i in data[2:3]:
                widths = [len(col)+1 for col in i.split()]

        for n in range(0, len(widths)):
            if n == len(widths)-1:
                colspecs.append((initial, -1))
            else:
                colspecs.append((initial, initial + widths[n]))
            initial +=  widths[n]

        # notice between bias and exp for best 1 and domain theres an extra space
        # deal with that here
        eleventh_column_start, eleventh_column_end = colspecs[10]
        eleventh_column_end += 1  # Increase the end index by 1 to widen the column
        colspecs[10] = (eleventh_column_start, eleventh_column_end)

        # Shift the subsequent columns
        for i in range(11, len(colspecs) - 1):
            start, end = colspecs[i]
            colspecs[i] = (start + 1, end + 1)

        # and only shift the start of the description
        start, end = colspecs[-1]
        colspecs[-1] = (start + 1, end)

        # finally apply the colspecs to the .tblout df
        hit_df = pd.read_fwf(hit_path, skipfooter=10, skiprows=3, colspecs = colspecs, names=['target_name', 'accession', 'query_name', 'accession_2', 'E-value', \
                                                            'score', 'bias', 'E-value_1_domain', 'score_1_domain', 'bias_1_domain', 'exp_domain', 'reg_domain', \
                                                              'clu_domain', 'ov_domain', 'env_domain', 'dom_domain', 'rep_domain', 'inc_domain', 'description_of_target'])

        if len(hit_df.target_name) < 0:
            print(row.target_name)
        else:
            starts = []
            ends = []
            strand_list  = []
            
            for i, hit in hit_df.iterrows():
                curr = parsed_fasta[hit.query_name]
                starts.append(curr[0])
                ends.append(curr[1])
                
                # if hit.accession == 'PF17733.4':
                #     print(hit.query_name)
                #     print(strands[hit.query_name])
                strand_list.append(strands[hit.query_name])

            hit_df['genome_start'] = starts
            hit_df['genome_end'] = ends

            hit_df['strand'] = strand_list
            hit_df['genome'] = [row.target_name for i in range(len(hit_df))]


            filt_hit_df = hit_df[hit_df['E-value'] < 0.01]
            # best_rows = filt_hit_df.groupby('query_name').apply(lambda x: x.nsmallest(3, 'E-value'))[['genome', 'accession', 'target_name', 'query_name', 'E-value', 'description_of_target', 'genome_start', 'genome_end', 'strand']]
            best_rows = filt_hit_df.groupby('query_name').apply(lambda x: x.nsmallest(1, 'E-value'))[['genome', 'accession', 'target_name', 'query_name', 'E-value', 'description_of_target', 'genome_start', 'genome_end', 'strand']]

            best_rows['intronID'] = row.intronID
            all_hit_df = all_hit_df.append(best_rows, ignore_index=True)


        all_hit_df['source'] = '.tblout'
        all_dom_df['source'] = '.domtblout'
        # all_hit_df['description_of_target'] = all_hit_df['description_of_target'].str.split(n=1).str[1]

    return pd.concat([all_hit_df, all_dom_df])


def combine_genbanks(all_hit_df, pharok_directory, bakta_directory, outdir, actual_genomes):
    """
    Run this to combine bakta and pharokka output
    """
    
    grouped_hits = all_hit_df.groupby('genome')

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Iterate over grouped hits
    for genome, hits_group in grouped_hits:
        if genome in actual_genomes:
            # Load GenBank file for the current genome
            pharokka_genbank_file = os.path.join(pharok_directory, f"{genome}/{genome}.gbk")
            bakta_genbank_file = os.path.join(bakta_directory, f"{genome}/{genome}.gbff")

            pharokka_genbank_record = SeqIO.read(pharokka_genbank_file, 'genbank')
            bakta_genbank_record = SeqIO.read(bakta_genbank_file, 'genbank')

            # Append features from the bakta record to the output file
            for feature in bakta_genbank_record.features:
                # print(feature)
                pharokka_genbank_record.features.append(feature)
            try:            
                SeqIO.write(pharokka_genbank_record, f"{outdir}/{genome}_double_annotated.gbff", 'genbank')
            except:
                print(f"error in genome {genome}")
 
 

           
def combine_defense(g2_df, bakta_dir, defense_dir):
    all_defense_hit_dfs = []

    for genome in g2_df.target_name.unique():
        map_df = pd.read_csv(f"{bakta_dir}/{genome}/{genome}.tsv", sep="\t", skiprows=5)
        df = pd.read_csv(f"{defense_dir}/{genome}/{genome}_defense_finder_hmmer.tsv", sep = '\t')

        best_rows = df.groupby('hit_id').apply(lambda x: x.nsmallest(2, 'i_eval')).reset_index(drop=True)

        for index, row in best_rows.iterrows():
            gene_row = map_df[map_df['Locus Tag'] == row['hit_id']]

            # Get the start and stop coordinates
            start = gene_row['Start'].iloc[0]
            stop = gene_row['Stop'].iloc[0]
            strand = gene_row['Strand'].iloc[0]

            aa_start = row.hit_begin_match
            aa_end = row.hit_end_match

            reading_frame = 1  # Assuming 1-based index

            # Convert amino acid positions to nucleotide positions
            hit_start_nt = (aa_start - 1) * 3 + reading_frame
            hit_end_nt = aa_end * 3 + reading_frame

            # Adjust for negative strand
            if strand == '-':
                hit_start_nt, hit_end_nt = hit_end_nt, hit_start_nt

            if strand == '-':
                hit_start_nt_absolute = stop - hit_start_nt + 1
                hit_end_nt_absolute = stop - hit_end_nt + 1
            else:
                hit_start_nt_absolute = start + hit_start_nt - 1
                hit_end_nt_absolute = start + hit_end_nt - 1

            best_rows.loc[index, 'hit_start_nt'] = int(hit_start_nt_absolute)
            best_rows.loc[index, 'hit_end_nt'] = hit_end_nt_absolute
            best_rows.loc[index, 'strand'] = strand


        best_rows['genome'] = genome

        all_defense_hit_dfs.append(best_rows)

    all_defense_hits = pd.concat(all_defense_hit_dfs, ignore_index=True)

    return all_defense_hits