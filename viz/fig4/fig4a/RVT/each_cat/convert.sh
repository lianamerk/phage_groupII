#!/usr/bin/env bash

for file in *.aout
do
  base="${file%.aout}"
  
  # 1) Convert .aout to .faa
  esl-reformat fasta "$file" > "${base}.faa"
  
  # 2) Use awk to take the first 20 sequences
  #    This command increments 'n' each time we see '>', 
  #    and prints lines only while n <= 20.
  awk '/^>/{n++} n <= 20 {print}' "${base}.faa" > "${base}_subset20.faa"
  
  echo "Processed $file -> ${base}.faa -> ${base}_subset20.faa"
done
