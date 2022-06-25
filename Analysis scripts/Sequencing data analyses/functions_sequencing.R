# Functions for Ribo-Seq and RNA-Seq analyses
#--------------------------------------------------
# Create colData for DESeq2 (convert strain number to strain name)
ccdt <- function(count_tab) {
  cdt <- data.frame(full_name = names(count_tab),
                    strain = sub("_.*", "", names(count_tab)),
                    size = sub(".*_", "", names(count_tab)),
                    CHX = ifelse(test = grepl(pattern = ".*minus.*", x = names(count_tab)), 
                                 yes = "minus", no = "plus"))
  cdt$ribosomes <- sub("^T.*", "TotalRNA", cdt$strain)
  cdt$ribosomes <- sub("^[QR][0-9].*", "TotalRiboseq", cdt$ribosomes)
  cdt$ribosomes <- sub("^[QR]U[0-9].*", "IPRiboseq", cdt$ribosomes)
  
  cdt$size <- sub("^T.*", "allnt", cdt$size)
  
  cdt$strain <- sub(".*504[123].*", "WTplusEV", cdt$strain)
  cdt$strain <- sub(".*533.*", "Upf1_FLAG_upf2del", cdt$strain)
  cdt$strain <- sub(".*503[345].*", "DE572AA_FLAG", cdt$strain)
  cdt$strain <- sub(".*501[7|9].*", "Upf1_FLAG", cdt$strain)
  cdt$strain <- sub(".*5020.*", "Upf1_FLAG", cdt$strain)
  cdt$strain <- sub(".*1[69].*", "FLAG_Upf1", cdt$strain)
  cdt$strain <- sub(".*69.*", "WT", cdt$strain)
  cdt$strain <- sub(".*70.*", "WT", cdt$strain)
  cdt$strain <- sub(".*490[23].*", "OE3", cdt$strain)
  
  cdt$lib <- paste0(cdt$strain, "@", cdt$CHX, "@", cdt$ribosomes)
  cdt$strainCHX <- paste0(cdt$strain, "@", cdt$CHX)
  cdt$ribosomes_size <- paste0(cdt$ribosomes, "@", cdt$size)
  cdt$ribosomesCHX <- paste0(cdt$ribosomes, "@", cdt$CHX)
  return(cdt)
}
#--------------------------------------------------
# Add mean of normalized counts to DESeq2 result table
norm_count <- function(dds, res, numerator_col, denominator_col, padj_cutoff = 0.05) {
  if (!("DESeq2" %in% (.packages()))) library(DESeq2)
  dds <- estimateSizeFactors(dds[, c(numerator_col, denominator_col)])
  #sizeFactors(dds)
  nrep_num <- length(numerator_col)
  nrep_denom <- length(denominator_col)
  normalized_counts <- as.data.frame(counts(dds, normalized = TRUE))
  norm_mean <- data.frame(row = row.names(normalized_counts),
                          ave_numerator = rowMeans(normalized_counts[, 1:nrep_num]),
                          ave_denominator = rowMeans(normalized_counts[, (nrep_num + 1):(nrep_num + nrep_denom)]))
  norm_mean <- full_join(norm_mean, res, by = "row")
  norm_mean$sig <- ifelse(test = (norm_mean$padj < padj_cutoff & !is.na(norm_mean$padj)), yes = paste0("padj < ", padj_cutoff), no = "ns")
  return(norm_mean)
}
#--------------------------------------------------
# Prepare table for plotting: add significance, group, NMD status, etc.
prepare_plot <- function(toplot, padj_cutoff = 0.05) {
  toplot$NMD <- ifelse(test = toplot$row %in% nmd_up, yes = "NMD substrates", no = "Non-NMD substrates")
  toplot$sig <- ifelse(test = (!is.na(toplot$padj) & toplot$padj < padj_cutoff), yes = paste0("padj < ", padj_cutoff), no = "ns")
  toplot$change <- "Unchanged"
  toplot[which(toplot$sig != "ns"), ]$change <- ifelse(test = toplot[which(toplot$sig != "ns"), ]$log2FoldChange > 0, yes = "Enriched", no = "Depleted")
  return(toplot)
}
#--------------------------------------------------