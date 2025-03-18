How to recreate this figure:

1. Search the sep2024 inphared database with the terminase LSU: PF03237
hmmsearch -o PF03237_millard_sept.out -A PF03237_millard_sept.aout --tblout PF03237_millard_sept.tblout PF03237.hmm /path/to/millard/genomes/11Sep2024_vConTACT2_proteins.faa 

2. Take the aout from this and construct an hmm: 
hmmbuild PF03237_millard.hmm PF03237_millard_sept.aout

3. Use this to search the IMG/VR genomes that have group II introns, and have been annotated with pharokka:
hmmsearch -o PF03237_millard_imgvr.out -A PF03237_millard_imgvr.aout --tblout PF03237_millard_imgvr.tblout millard/PF03237_millard.hmm /path/to/IMGVR/genomes/pharokka_annotation/phanotate.faa
and also to search the Millard database again:
hmmsearch -o PF03237_millard_sept_2.out -A PF03237_millard_sept_2.aout --tblout PF03237_millard_sept_2.tblout millard/PF03237_millard.hmm /path/to/millard/genomes/11Sep2024_vConTACT2_proteins.faa  

4. Run through the process in millard/terl_millard_processing.ipynb. This will eventually give you a file, millard/Terminase_6N.faa, which you will use in the joint terminase tree.

5. Open millard_and_imgvr/terl_joint_processing.ipynb and follow through, ending with Terminase_6N_imgvr_cluster_rep_single_name.faa

6. Concat prior tree and the IMGVR tree (millard/Terminase_6N.faa and millard_and_imgvr/Terminase_6N_imgvr_cluster_rep_single_name.faa) > prior_tree_and_imgvr_reps.faa

7. Run tree_together.sh