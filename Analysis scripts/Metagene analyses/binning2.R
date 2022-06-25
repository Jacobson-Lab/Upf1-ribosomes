# Function to determine footprint abundance in bins of normalized coding region
  # For individual gene
    # dat: a data sample in reads_psite_list 
    # bins: number of bins to split coding region

binning2 <- function(dat, bins = 4) {
  dat$gene_length <- dat$cds_stop - dat$cds_start + 1
  cds_sub <- dat[psite_region == "cds", c("transcript", "gene_length", "psite_from_start")]
  cds_sub$bin <- floor((cds_sub$psite_from_start)/(cds_sub$gene_length/bins))
  count_tab <- cds_sub[, list(N = .N), by = list(bin, transcript, gene_length)]
  count_tab <- reshape2::dcast(data = count_tab, transcript + gene_length ~ bin, value.var = "N")
  count_tab[is.na(count_tab)] <- 0
  count_tab$bin_total <- rowSums(count_tab[, -c(1:2)])
  count_tab2 <- cbind(count_tab[, c(1:2, length(count_tab))], count_tab[, -c(1:2, length(count_tab))] / count_tab$bin_total)
  count_tab2$occupancy <- count_tab2$bin_total / count_tab2$gene_length
  count_tab2 <- reshape::melt(count_tab2, id.vars = c("transcript", "gene_length", "bin_total", "occupancy"))
  return(count_tab2)
}