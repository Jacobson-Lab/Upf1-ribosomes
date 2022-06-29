# ---------------------------------------
# Calculate fraction of footprints P-site in each region and reading frame
  # Use only footprint lengths 20-23nt, 27-32nt, 37-43nt, with corrected/verified P-site offsets
# ---------------------------------------

# Fraction of footprints in each reading frame, overall or by region --> Data for Figure S8 C-D
psite_fraction <- function(dat) {
  df <- data.frame(table(psite_region = dat$psite_region, Frame = dat$frame))
  df$fraction_overall <- df$Freq/sum(df$Freq)
  reg_freq <- data.frame(table(psite_region = dat$psite_region))
  df <- dplyr::left_join(df, reg_freq, by = "psite_region", suffix = c("", "_region_sum"))
  df$fraction_by_region <- df$Freq/df$Freq_region_sum
  return(df)
}

pf_list <- lapply(reads_psite_list, psite_fraction)
pf <- bind_rows(pf_list, .id = "sample")

# Fraction of footprints in each region, regardless of reading frame --> Data for Figure S8 A-B
rf_list <- lapply(pf_list, function(x) {
  x <- unique(x[, c("psite_region", "Freq_region_sum")])
  x$fraction <- x$Freq_region_sum/sum(x$Freq_region_sum)
  return(x)
})
rf <- bind_rows(rf_list, .id = "sample")
