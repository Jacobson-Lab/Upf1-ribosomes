#!/bin/bash
# Automate automate_codon_window_count_by_codonpos_from_rpl.sh
#	Argument 1: path to input files
# 	Argument 2: path to output files
#   Argument 3: path to txt file containing list of codons to perform analysis

files=$(find $1 -name "*.txt")
for f in $files
do  fname=$(basename "$f")
	fname2="${fname%.txt}"
	path_out="$2/$fname2"
	bsub -q short -W 3:00 -R "rusage[mem=3000]" -R "span[hosts=1]" -n 1 ./automate_codon_window_count_by_codonpos_from_rpl.sh $f $path_out $3
done