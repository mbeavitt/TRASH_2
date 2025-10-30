sequence_window_score <- function(fasta_sequence, window_size, kmer = 10) {
  # TODO check if changing these settings below can make the script work better,
  # although these were optimised
  if ((window_size / 2) <= kmer) stop("sequence_window_score: window size is too small")

  sequence_full_length <- length(fasta_sequence)

  result <- genomic_bins(start = 1, end = sequence_full_length, bin_size = window_size)
  starts <- result$starts
  ends <- result$ends

  # Process all windows in one go
  scores <- NULL
  cat("Sequence full length: ", round(sequence_full_length / 1000000, 3), " Mbp\n", sep = "")
  cat("Processing ", length(starts), " windows...\n", sep = "")

  for (i in seq_along(starts)) {
    result <- seq_win_score_int(1, window_size, kmer, fasta_sequence[starts[i] : ends[i]])
    scores <- c(scores, result)
  }

  cat("Done!\n")

  if ((sum(is.na(scores))>0) || (length(scores) == 0)) {
    stop("sequence_window_score did not produce valid result")
  }
  return(scores)
}
