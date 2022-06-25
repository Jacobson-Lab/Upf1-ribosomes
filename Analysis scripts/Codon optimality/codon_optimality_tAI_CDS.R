#----------------------------------------------
# Calculate tAI (tRNA adaptation index) of CDS
# ---------------------------------------------
library(Biostrings)
library(seqinr)
library(dplyr)

## Get CDS sequence for each mRNA entry
# Read in sequence data
fasta <- readDNAStringSet("https://raw.githubusercontent.com/Jacobson-Lab/yeast_transcriptome_v5/main/gffread_transcripts_from_v5_transcriptome_13AUG20218.fa")
df <- data.frame(transcript_full = names(fasta), sequence = paste(fasta))
df$transcript_id <- sub(" gene.*", "", df$transcript_full)

# Add sequence lengths data and extract only CDS sequence
transcriptome <- read.table("../Sequencing data analyses/transcriptome_v5_mRNA_region_length.txt", header = TRUE)
df <- left_join(df, transcriptome[, c("transcript_id", "length", "l_utr5", "l_cds", "l_utr3")], by = "transcript_id")
df[is.na(df$l_cds), ]$l_cds <- df[is.na(df$l_cds), ]$length
df[is.na(df$l_utr5), ]$l_utr5 <- 0
df[is.na(df$l_utr3), ]$l_utr3 <- 0
seq_cds <- substr(df$sequence, start = df$l_utr5 + 1, stop = df$l_utr5 + df$l_cds)

## Codon optimality score (a.k.a "w") for each codon
# Derived based on criteria developed by dos Reis et al., 2004 (doi: 10.1093/nar/gkh834).
# Obtained the calculated values from Tuller et al., 2010 (doi: 10.1073/pnas.0909910107)
scw <- read.table("codon_optimality_w.txt", header = TRUE)

## Functions
# Function to tabulate codon frequency, excluding stop codons
codon_c <- function(s) { 
  x <- as.data.frame(uco(s2c(s), frame = 0, index = "eff"))
  x <- x[which(!(x$Var1 %in% c("taa", "tag", "tga"))), ]
  rownames(x) <- as.character(toupper(unlist(x[, "Var1"])))
  x <- x[match(scw$codons, (rownames(x))), ]
  x <- matrix(t(x[, 2]), ncol = 61)
  return(x)
}
# Function to automate codon_c over a list of sequences
codonMatrix <- function(data) { 
  mat <- matrix(0L, nrow = length(data), ncol = 61)
  for (i in 1:length(data)) {
    mat[i, ] <- codon_c(data[i])
  }
  return(mat)
}
# Function to calculate tAI. From https://github.com/mariodosreis/tai
get.tai <- function(x, w) {
  w = log(w)              # calculate log of w
  n = apply(x,1,'*',w)    # multiply each row of by the weights
  n = t(n)                # transpose
  n = apply(n,1,sum)      # sum rows
  L = apply(x,1,sum)      # get each ORF length
  tAI = exp(n/L)          # get tai
  return(tAI)
}

## Get tAI
mat_cds <- codonMatrix(seq_cds)
tAI_cds <- get.tai(x = mat_cds, w = scw$w)
opt_df <- data.frame(transcript_id = transcriptome[, "transcript_id"], tAI_cds = tAI_cds)
write.table(opt_df, file = "transcriptome_v5_codon_optimality_tAI_cds.txt", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
