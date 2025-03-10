How to recreate this figure:

1. Text files with each of the Toro categories are in the directory each_cat (for retron, gii, g2l, crispr, dgr, abi). First pull those header prefixes from the bigger 9141_toro.faa using the following commands:
python pull.py 9141_toro.faa retron.txt retron.faa
python pull.py 9141_toro.faa DGR.txt DGR.faa
python pull.py 9141_toro.faa abi.txt abi.faa
python pull.py 9141_toro.faa GII.txt GII.faa
python pull.py 9141_toro.faa crispr.txt crispr.faa
python pull.py 9141_toro.faa g2l.txt g2l.faa

2. Then, hmmsearch using each fasta file and the Pfam RVT model:
hmmsearch --tblout RVT_retron_hmmscan.tblout -o RVT_retron_hmmscan.out -A RVT_retron_hmmscan.aout ./PF00078.hmm retron.faa
hmmsearch --tblout RVT_DGR_hmmscan.tblout -o RVT_DGR_hmmscan.out -A RVT_DGR_hmmscan.aout ../PF00078.hmm DGR.faa
hmmsearch --tblout RVT_abi_hmmscan.tblout -o RVT_abi_hmmscan.out -A RVT_abi_hmmscan.aout ../PF00078.hmm abi.faa
hmmsearch --tblout RVT_gii_hmmscan.tblout -o RVT_gii_hmmscan.out -A RVT_gii_hmmscan.aout ../PF00078.hmm GII.faa
hmmsearch --tblout RVT_crispr_hmmscan.tblout -o RVT_crispr_hmmscan.out -A RVT_crispr_hmmscan.aout ../PF00078.hmm crispr.faa
hmmsearch --tblout RVT_g2l_hmmscan.tblout -o RVT_g2l_hmmscan.out -A RVT_g2l_hmmscan.aout ../PF00078.hmm g2l.faa

3. Convert each .out to a .faa and take a subset of 20 of those so the tree isn't too full. The script to do this is called "convert.sh" so just run the following command in the each_cat directory.
bash convert.sh

4. Combine them into one big .faa: cat *subset20.faa > toro_evensubset.faa.

Now, to collect the sequences from the millard introns, (still in the each_cat directory), run:
sbatch -t 6-00:00 -p eddy -J hmmsearch -c 64 --mem-per-cpu=6G -N 1 -o %x.out --wrap="hmmsearch --tblout RVT_millard_hmmscan.tblout -o RVT_millard_hmmscan.out -A RVT_millard_hmmscan.aout ./PF00078.hmm /n/eddy_lab/users/lmerk/millard_sept/inphared_11Sep2024/11Sep2024_vConTACT2_proteins.faa"

5. From the Millard search, grab just the ones that are in group II introns: python pull.py millard_RVT_hmmsearch.faa  millard_introns.txt millard_intron.faa

6. Combine the two: cat toro_evensubset.faa millard_intron.faa > both.faa

7. Submit a job to make the tree (in an environment that has MAFFT, mine is in 'defensefinder' env):
sbatch tree.sh

8. View it in iTol, and you can use the metadata text file (final_colorstrip.txt).

