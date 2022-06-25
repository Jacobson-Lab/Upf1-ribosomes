# Function to tabulate read length and calculate fractions of read length in a ribo-seq library
    # dt: an input data table in reads_list or reads_psite_list

rl_dist <- function(dt) {
  rl <- data.table(table(dt$length))
  names(rl) <- c("length", "count")
  rl$fraction <- rl$count/sum(rl$count)
  return(rl)
}