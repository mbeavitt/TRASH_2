genomic_bins <- function(start = 1, end = 0, bin_size = 0) {
  if (end < start) stop("genomic bins: End smaller than start will not work too well...")

  # Handle case where we only need one bin
  if (end <= bin_size) {
    return(list(starts = start, ends = end))
  }

  # Calculate start positions
  start_positions <- seq(start, (end - bin_size), bin_size)

  # Remove last window if it would be less than half a bin size
  if ((end - start_positions[length(start_positions)]) < (bin_size / 2)) {
    start_positions <- start_positions[-length(start_positions)]
  }

  # Calculate end positions
  if (length(start_positions) == 1) {
    end_positions <- end
  } else {
    end_positions <- c((start_positions[-1] - 1), end) + bin_size
  }

  # Cap ends at sequence length
  end_positions[end_positions > end] <- end

  return(list(starts = start_positions, ends = end_positions))
}

genomic_bins_starts <- function(start = 1, end = 0, bin_number = 0, bin_size = 0) {
  if (bin_number > 0 && bin_size > 0) stop("genomic bins starts: Use either bin number or bin size")
  if (bin_number == 0 && bin_size == 0) stop("genomic bins starts: Use either bin number or bin size")
  if (end < start) stop("genomic bins starts: End smaller than start will not work too well...")
  if (bin_size >= (end - start)) return(start)

  if (bin_number > 0) {
    seq_per_bin <- (end - start + 1) %/% bin_number
    remaining_seq <- (end - start + 1) %% bin_number
    bin_sizes <- rep(seq_per_bin, bin_number)
    if (remaining_seq > 0) {
      add_remaining_here <- sample(1:bin_number, remaining_seq)
      bin_sizes[add_remaining_here] <- bin_sizes[add_remaining_here] + 1
    }
    start_positions <- bin_sizes
    for (i in 2 : length(start_positions)) {
      start_positions[i] <- start_positions[i] + start_positions[i - 1]
    }
    start_positions <- start_positions - start_positions[1] + start
    remove(seq_per_bin, remaining_seq, bin_sizes)
    return(start_positions)
  }
  if (bin_size > 0) {
    if ((end - start) < bin_size) return(start)
    start_positions <- seq(start, (end - bin_size), bin_size)
    if ((end - start_positions[length(start_positions)]) < (bin_size / 2)) { #if future last win length is less than half a bin size
      start_positions <- start_positions[-length(start_positions)] # remove the last one
    }
    return(start_positions)
  }
  return(NA)
}
