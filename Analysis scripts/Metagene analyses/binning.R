# Function to determine footprint abundance in bins of normalized coding region
  # Meta-gene
    # dat: a data sample in reads_psite_list 
    # bins: number of bins to split coding region
    # transcripts: a vector containing name of genes of interest; If NULL, all genes are considered.

binning <- function(dat, bins = 100, transcripts = NULL) {
  dat$gene_length <- dat$cds_stop - dat$cds_start + 1
  if (length(transcripts) == 0) {
    cds_sub <- dat[psite_region == "cds", c("transcript", "gene_length", "psite_from_start")]
  } else {
    cds_sub <- dat[psite_region == "cds" & transcript %in% transcripts, c("transcript", "gene_length", "psite_from_start")]
  }
  lib_size <- nrow(cds_sub)
  cds_sub$bin <- floor((cds_sub$psite_from_start+1)/(cds_sub$gene_length/bins))
  count_tab <- cds_sub[, list(N = .N), by = list(bin)]
  count_tab$fraction <- count_tab$N/lib_size
  return(count_tab)
}