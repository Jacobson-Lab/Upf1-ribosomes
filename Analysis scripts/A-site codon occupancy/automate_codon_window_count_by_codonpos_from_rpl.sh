#!/bin/bash
# Automate codon_window_count_by_codonpos_from_rpl.py
#	Argument 1: input file
# 	Argument 2: base name for output file
#   Argument 3: path to txt file containing list of codons to perform analysis

module load python3/3.5.0
module load R/3.5.0

filename=$(basename "$1")
filename2="${filename%.txt}"

mkdir temp_$filename2
codon_list=$(cat $3)
for codon in $codon_list
do  echo $codon
	python3 codon_window_count_by_codonpos_from_rpl.py $1 ./temp_$filename2/temp A $codon 0 30
	Rscript --vanilla calc_REV_window_codonpos_v2.R ./temp_$filename2/temp_$codon.txt 30 0.1
	rm ./temp_$filename2/temp_$codon.txt
done

files=$(find ./temp_$filename2/ -name "temp_nt*.txt")
for k in $files
do  full_path_out="$2_bynt.txt"
	echo $k
	cat $k >> $full_path_out
	rm $k
done

files=$(find ./temp_$filename2/ -name "temp_codon*.txt")
for k in $files
do  full_path_out="$2_bycodon.txt"
	echo $k
	cat $k >> $full_path_out
	rm $k
done

rmdir temp_$filename2