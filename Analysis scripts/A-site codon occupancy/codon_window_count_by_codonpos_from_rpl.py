# Count reads within a window around a specified codon for each gene, calculate relative enrichment
#	Argument 1: input file
# 	Argument 2: base name for output file
#	Argument 3: E, P, or A site
#	Argument 4: codon of interest (individual or pair or any length)
#	Argument 5: reading frame (0, 1, or 2)
#	Argument 6: window around the codon of interest (e.g. 50 for 50 codons around the codon of interest)

import sys
import pickle

# Load pickle file containing sequence of mRNAs in the transcriptome
tp = open("transcriptome_v5.pickle", "rb")
seq_dict = pickle.load(tp)
tp.close()

# Function to identify positions of codon of interest in a sequence. Return the position of the first nucleotide of that codon.
def codon_pos(codon, gene, start, end, frame):
	sequence = seq_dict[gene][start:end]
	#print(sequence)
	#print(len(sequence)/3)
	positions = []
	l = len(codon)
	if codon in sequence:
		for i in range(frame, len(sequence), 3): # Loop every 3 letters, starting at specified frame (0, 1, or 2):
			x = sequence[i:i+l]
			if x == codon:
				positions.append(i)
	return positions

# Function to count reads within a window for a codon
def count_in_window(cod, infile, frame, site, window):
	# Initiating results dictionaries
	pos_dict = {} 	 # Dictionary storing positions of the codon of interest within each gene [key]
	count_dict = {}  # Dictionary storing read count within the window for each position [key]
	# infile's column: 0-transcript 1-size 2-psite 3-cds_start 4-cds_stop 5-psite_from_start 6-psite_from_stop 7-psite_region 8-frame 9-count
	for line in open(infile):
		line_split = line.strip().split('\t')
		transcript = line_split[0]
		start = int(line_split[3])-1 # minus 1 because python starts at 0 while R starts at 1
		end = int(line_split[4])+3 # add 3 to include stop codon
		if transcript not in pos_dict:
			pos_dict[transcript] = codon_pos(cod, transcript, start, end, frame)
		if site == 'E':
			pos = int(line_split[5])-3
		elif site == 'P':
			pos = int(line_split[5])
		elif site == 'A':
			pos = int(line_split[5])+3
		else:
			print("Invalid site.")
			break
		for x in pos_dict[transcript]:
			if x > 0: # Do not consider start codon
				dist = pos-x # Distance of read to position of codon of interest
				if abs(dist) <= window: # If distance is within window, add read count to that distance
					dpos = dist + window # convert distance to position in list
					c = int(line_split[9])
					if str(line_split[1]) == 'M':
						s = 'M'
					else:
						s = 'Z' # L or S
					cpos = transcript + '_' + str(x) + '_' + s
					if cpos not in count_dict:
						count_dict[cpos] = [0] * ((window*2)+1) # Initialize list of 0's 
					count_dict[cpos][dpos] += c
	return count_dict

# Initiating input variables
infile = str(sys.argv[1])
outfile = str(sys.argv[2])
site = str(sys.argv[3])
cod = str(sys.argv[4]) 
f = int(sys.argv[5])
window = int(sys.argv[6]) * 3 # convert to nt

res = count_in_window(cod, infile, f, site, window)
out = open(outfile + "_" + str(cod) + ".txt", "w")
for cpos, dpos in res.items():
	out.write(str(cpos) + '\t' + '\t'.join(str(c) for c in dpos) + '\n')
