genomic_bins_starts <- function(start = 1, end = 0, bin_size = 0) {
  if (end < start) stop("genomic bins starts: End smaller than start will not work too well...")
  if ((end - start) <= bin_size) return(start)

  start_positions <- seq(start, (end - bin_size), bin_size)
  if ((end - start_positions[length(start_positions)]) < (bin_size / 2)) { #if future last win length is less than half a bin size
    start_positions <- start_positions[-length(start_positions)] # remove the last one
  }
  return(start_positions)
}
