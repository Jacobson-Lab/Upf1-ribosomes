#----------------------------------------------
# Prepare a data table in riboWaltz's reads_psite_list object for A-site occupancy analysis
  # Add footprint size information
  # Collapse reads with the same psite location into a single entry and add the "count" column
# ---------------------------------------------

# dt: data table in reads_psite_list
# samp: name of sample
# CHX: Whether sample was treated with CHX. -CHX = "minus" (default), +CHX = "plus". Different footprint lengths are retained for each case.
# folder: path to the folder where the output file will be written to. (The folder must be created first)

combine_identical <- function(dt, samp, CHX = "minus", folder = "./rpl/") {
  require(data.table)
  if (CHX == "minus") {
    dt <- dt[length %in% 20:32, ]
    dt$size <- ifelse(test = dt$length > 25, yes = "M", no = "S")
  } else if (CHX == "plus") {
    dt <- dt[length %in% 27:43, ]
    dt$size <- ifelse(test = dt$length < 33, yes = "M", no = "L")
  } else {
    message("Invalid CHX option. Specify 'plus' or 'minus'. Default is 'minus'")
  }
  dt <- dt[, c("transcript", "size", "psite", "cds_start", "cds_stop", "psite_from_start", "psite_from_stop", "psite_region", "frame")]
  xx <- unique(dt[ , count := .N, by = list(transcript, size, psite)])
  write.table(xx, file = paste0(folder, samp, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

# Single table:
#combine_identical(dt = reads_psite_list[["WT"]], samp = "WT", CHX = "plus")

# Multiple tables from a list of data tables
sapply(names(reads_psite_list), function(x) {combine_identical(dt = reads_psite_list[[x]], samp = x, CHX = "plus")})
