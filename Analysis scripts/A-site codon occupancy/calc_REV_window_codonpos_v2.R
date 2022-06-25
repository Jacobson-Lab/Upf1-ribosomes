# Calculate relative occupancy in a window
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
window <- as.numeric(args[2]) # Window in codons
ao_cutoff <- as.numeric(args[3])
codon <- sub(".txt", "", infile)
codon <- sub(".*_", "", codon)
path <- sub(paste0("/temp_", codon, ".txt"), "", infile)

df <- read.table(infile, header = FALSE)
message("Read table")
colnames(df) <- c("pos", seq(-(window*3), window*3))

# Calculate number of reads per codon in an mRNA to use for filtering minimum occupancy
df$total <- rowSums(df[, -c(1)])
df$ave_occupancy <- df$total / ((window*2) + 1)
df$ave_nt_occupancy <- df$total / (((window*2) + 1)*3)
df_filter <- df[which(df$ave_occupancy > ao_cutoff), ] # only keep mRNAs passing minimum cutoff. E.g. 0.2 would mean 1 read per 5 codons
message("Finish filtering\n", "Original: ", nrow(df), ", Filtered: ", nrow(df_filter))
rm(df)

# Prepare data
notpos <- c("pos", "total", "ave_occupancy", "ave_nt_occupancy")
df_cols <- df_filter[, notpos]
df_pos <- df_filter[, which(!(colnames(df_filter) %in% notpos))]
rm(df_filter)

# Mean codon occupancy: nt resolution ----------------------------------------------------------------------------------------
message("Processing nt resolution")
df_norm <- apply(df_pos, 2, "/", df_cols$ave_nt_occupancy)
indexM <- grepl(pattern = "_M$", x = df_cols$pos)
if (sum(indexM) >= 2) {
  df_M <- data.frame(mean_occ = colMeans(df_norm[indexM, ]))
  df_M$Distance <- as.numeric(rownames(df_M))
} else {
  df_M <- data.frame(mean_occ = NA, Distance = NA)
}
indexZ <- grepl(pattern = "_Z$", x = df_cols$pos)
if (sum(indexZ) >= 2){
  df_Z <- data.frame(mean_occ = colMeans(df_norm[indexZ, ]))
  df_Z$Distance <- as.numeric(rownames(df_Z))
} else {
  df_Z <- data.frame(mean_occ = NA, Distance = NA)
}
df_nt <- dplyr::bind_rows(list(M = df_M, Z = df_Z), .id = "size")
df_nt <- df_nt[!is.na(df_nt$Distance), ]
df_nt$codon <- codon
write.table(df_nt, file = paste0(path, "/temp_nt_", codon, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
rm(df_norm, df_M, df_Z, df_nt)

# Mean codon occupancy: codon resolution -------------------------------------------------------------------------------------
message("Processing codon resolution")
dff <- df_pos
colnames(dff) <- floor(as.numeric(colnames(dff)) / 3) # Rename distance from nt to codon resolution
dff2 <- as.data.frame(do.call(cbind, lapply(split(seq_len(ncol(dff)), colnames(dff)), function(x) rowSums(dff[x])))) # Combine counts from the same codon position
df_norm <- apply(dff2, 2, "/", df_cols$ave_occupancy)
indexM <- grepl(pattern = "_M$", x = df_cols$pos)
if (sum(indexM) >= 2) {
  df_M <- data.frame(mean_occ = colMeans(df_norm[indexM, ]))
  df_M$Distance <- as.numeric(rownames(df_M))
} else {
  df_M <- data.frame(mean_occ = NA, Distance = NA)
}
indexZ <- grepl(pattern = "_Z$", x = df_cols$pos)
if (sum(indexZ) >= 2){
  df_Z <- data.frame(mean_occ = colMeans(df_norm[indexZ, ]))
  df_Z$Distance <- as.numeric(rownames(df_Z))
} else {
  df_Z <- data.frame(mean_occ = NA, Distance = NA)
}
df_cd <- dplyr::bind_rows(list(M = df_M, Z = df_Z), .id = "size")
df_cd <- df_cd[!is.na(df_cd$Distance), ]
df_cd$codon <- codon
write.table(df_cd, file = paste0(path, "/temp_codon_", codon, ".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
rm(df_norm, df_M, df_Z, df_cd)


