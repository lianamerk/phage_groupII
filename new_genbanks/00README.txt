This directory implements Infernal and updates the genbank files with the individual Infernal hits (in genbank_with_intron_hits). Then, once a curated version of the intron dataframe is obtained (with full length start/ends, saved in g2_df.csv), add_full_length_intron.py can be run to add the whole intron annotation to the genbank file (in genbank_with_full_intron directory).

1. Run the bash script to submit Infernal jobs `infernal.sh`
2. Run `python main.py`
3. Once you have curated the start/ends of the introns into g2_df.csv, you can run `python add_full_length_intron.py`
