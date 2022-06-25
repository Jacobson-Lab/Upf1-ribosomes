# Functions for mass spec data analyses
#--------------------------------------------------
# Function to consolidate iBAQ, Protein Identification Probability, and Total Spectrum Count information (exported from Scaffold)
#   Criteria: 1) Protein identification probability >= 99% --> pass filter
#             2) Total spectrum count >= 2 --> pass filter
consolidate <- function(filepath) {
  require(readxl)
  dat <- read_xlsx(path = filepath, col_names = FALSE, skip = 2)
  cols <- read_xlsx(path = filepath, col_names = FALSE, n_max = 2)
  # Split information
  ibaq <- dat[, grepl(pattern = "iBAQ", x = cols[1, ])]
  ibaq[ibaq == -1] <- 0 # Replace -1 with 0
  pip <- dat[, grepl(pattern = "Probability", x = cols[1, ])]
  tsc <- dat[, grepl(pattern = "Spectrum", x = cols[1, ])]
  # Filter
  pip_pass <- pip >= 0.99
  tsc_pass <- tsc >= 2
  ibaq_pass <- ibaq * pip_pass # Make iBAQ of pip_pass == FALSE equals 0
  ibaq_pass <- ibaq_pass * tsc_pass # Make iBAQ of tsc_pass == FALSE equals 0
  ibaq_pass <- cbind(dat[, grepl(pattern = "Alternate ID", x = cols[2, ])], ibaq_pass)
  colnames(ibaq_pass) <- cols[2, colnames(ibaq_pass)]
  return(ibaq_pass)
}
#--------------------------------------------------
# Function to remove rows with zero values among replicates and average replicates
ave_nonzero <- function(dat, cols, minrep = 2) { # minrep = minimum of replicates with non-zero iBAQ
  df <- dat[, c(1, cols)]
  keep <- rowSums(df[, -1] > 0) >= minrep
  df <- df[keep, ]
  df2 <- df
  df2[, "iBAQ"] <- rowMeans(df2[, -1])
  df2 <- left_join(dat[, "Alternate ID", drop = FALSE], df2[, c("Alternate ID", "iBAQ")], by = "Alternate ID")
  df_list <- list(separate = df, average = df2)
  return(df_list)
}
#--------------------------------------------------
# Wrapper function for ave_nonzero, do a pair of data (IP and Total) and combine them into one data.frame
ave_nonzero_combine <- function(dat, Total_col, IP_col, minrep = 2) {
  Total_list <- ave_nonzero(dat = dat, cols = Total_col, minrep = minrep)
  IP_list <- ave_nonzero(dat = dat, cols = IP_col, minrep = minrep)
  # Separate df, for limma
  df <- full_join(Total_list$separate, IP_list$separate, by = "Alternate ID")
  colnames(df) <- c("Protein", paste0("Total_", 1:(length(Total_list$separate)-1)), paste0("IP_", 1:(length(IP_list$separate)-1)))
  # Average df, for scatter plot
  df2 <- full_join(Total_list$average, IP_list$average, by = "Alternate ID")
  colnames(df2) <- c("Protein", "Total", "IP")
  df_list <- list(separate = df, average = df2)
  return(df_list)
}
#--------------------------------------------------
# Function to impute missing values
#   Adapted from https://github.com/wasimaftab/LIMMA-pipeline-proteomics
#   Mimic: Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
#   data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3) 
#   https://bioconductor.org/packages/devel/bioc/vignettes/DEP/inst/doc/DEP.html#workflow-functions-for-the-entire-analysis
#   https://www.nature.com/articles/nprot.2017.147
#   Imputation from this paper is based on MSnbase R package, but this "man" option is no longer available
impute_na <- function(dat) { # Make sure input is normally distributed (e.g. log-transform data before using this function)
  data_limma <- as.matrix(dat)
  data_limma[is.infinite(data_limma)] <- NA
  nan_idx <- which(is.na(data_limma))
  fit <- fitdistr(c(na.exclude(data_limma)), "normal")
  mu <- as.double(fit$estimate[1])
  sigma <- as.double(fit$estimate[2])
  sigma_cutoff <- 6
  new_width_cutoff <- 0.3
  downshift <- 1.8
  width <- sigma_cutoff * sigma
  new_width <- width * new_width_cutoff
  new_sigma <- new_width / sigma_cutoff
  new_mean <- mu - downshift * sigma
  imputed_vals_my = rnorm(length(nan_idx), new_mean, new_sigma)
  data_limma[nan_idx] <- imputed_vals_my
  return(data_limma)
}
#--------------------------------------------------
# Function to run imputation for missing values and fit model for differential expression between IP and Total
eb_fit <- function(ibaq_pass) {
  # Identify exclusive proteins in Total
  n_IP_rep <- ncol(ibaq_pass[, grepl(pattern = "IP", x = colnames(ibaq_pass))])
  IP_NA <- rowSums(is.na(ibaq_pass[, grepl(pattern = "IP", x = colnames(ibaq_pass))])) == n_IP_rep
  Total_exclusive <- ibaq_pass[IP_NA, ]
  # Identify exclusive proteins in IP
  n_Total_rep <- ncol(ibaq_pass[, grepl(pattern = "Total", x = colnames(ibaq_pass))])
  Total_NA <- rowSums(is.na(ibaq_pass[, grepl(pattern = "Total", x = colnames(ibaq_pass))])) == n_Total_rep
  IP_exclusive <- ibaq_pass[Total_NA, ]
  # Filter out exclusive proteins
  rownames(ibaq_pass) <- ibaq_pass$Protein
  ibaq_filter <- ibaq_pass[!(ibaq_pass$Protein %in% c(Total_exclusive$Protein, IP_exclusive$Protein)), ]
  ibaq_log <- log2(as.matrix(ibaq_filter[, -1]))
  # Impute missing values
  ibaq_log <- impute_na(ibaq_log)
  # Set up design matrix
  xdesign <- data.frame(fullname = colnames(ibaq_log))
  xdesign$ribo <- ifelse(test = grepl(pattern = "Total", x = xdesign$fullname), yes = "Total", no = "IP")
  group <- factor(xdesign$ribo)
  group <- factor(group, levels = c("Total", "IP"))
  design_matrix <- model.matrix(~0+group)
  # Set up contrast matrix
  cont <- makeContrasts(IPvsTotal = groupIP-groupTotal, levels = colnames(design_matrix))
  # Run limma
  fit <- lmFit(ibaq_log, design_matrix)
  fit <- contrasts.fit(fit, contrasts = cont)
  fit.eb <- eBayes(fit)
  res <- topTable(fit.eb, number = Inf, adjust.method = "BH")
  res$Protein <- rownames(res)
  # Return results
  res_list <- list(IP_exclusive = IP_exclusive, Total_exclusive = Total_exclusive, filtered = ibaq_filter, eb_fit_res = res)
  return(res_list)
}
#--------------------------------------------------
# Function to add protein group labels to data table
label_protein <- function(res_df) {
  # Load/define list of proteins of interest
  # Ribosomal proteins
  rpsl <- scan("./protein_group_list/ribosomal_protein_genes.txt", character())
  # initiation factors
  ifg <- scan("./protein_group_list/translation_initiation_factor_genes.txt", character())
  # elongation factors
  efg <- scan("./protein_group_list/translation_elongation_factor_genes.txt", character())
  # termination factors
  tfg <- c("SUP45", "SUP35")
  # decay factors
  dfg <- scan("./protein_group_list/decay_factor_genes.txt", character())
  # protein folding chaperone
  pcg <- scan("./protein_group_list/protein_folding_chaperone_genes.txt", character())
  
  # For point color
  res_df$Protein2 <- "Other"
  tryCatch({res_df[which(res_df$Protein == "NAM7"), ]$Protein2 <- "Upf1"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein == "NMD2"), ]$Protein2 <- "Upf2"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein == "UPF3"), ]$Protein2 <- "Upf3"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein %in% rpsl), ]$Protein2 <- "Ribosomal protein"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein %in% ifg), ]$Protein2 <- "Translation initiation factor"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein %in% efg), ]$Protein2 <- "Translation elongation factor"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein %in% tfg), ]$Protein2 <- "Translation termination factor"}, error = function(e){}) 
  tryCatch({res_df[which(res_df$Protein %in% dfg), ]$Protein2 <- "mRNA decay factor"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein %in% pcg), ]$Protein2 <- "Protein folding chaperone"}, error = function(e){})
    
  # For label
  tryCatch({res_df$repel <- ""}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein == "NAM7"), ]$repel <- "Upf1"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein == "NMD2"), ]$repel <- "Upf2"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein == "UPF3"), ]$repel <- "Upf3"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein == "SUP45"), ]$repel <- "eRF1"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein == "SUP35"), ]$repel <- "eRF3"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein == "DCP1"), ]$repel <- "Dcp1"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein == "DCP2"), ]$repel <- "Dcp2"}, error = function(e){})
  tryCatch({res_df[which(res_df$Protein == "EDC3"), ]$repel <- "Edc3"}, error = function(e){})

  return(res_df)
}
#--------------------------------------------------
# Function to combine data and export them into csv format
export_data <- function(data_list, sample_name) {
  ibaqs <- bind_rows(data_list[1:3], .id = "Detection")
  ibaqs$Detection <- recode(ibaqs$Detection, filtered = "Detected_in_both")
  write.csv(ibaqs, file = paste0("./export_data/", sample_name, "_filtered_iBAQ_values.csv"), row.names = FALSE)
  if (length(data_list) == 4) {
    eb <- data_list[[4]]
    write.csv(eb, file = paste0("./export_data/", sample_name, "_limma_abundance_changes.csv"), row.names = FALSE)
  }
}
#--------------------------------------------------