sequence_window_score <- function(fasta_sequence, window_size, kmer = 10, output_dir = ".") {
  # TODO check if changing these settings below can make the script work better,
  # although these were optimised
  fraction_p <- 0.5
  if ((window_size / 2) <= kmer) stop("sequence_window_score: window size is too small")

  sequence_full_length <- length(fasta_sequence)

  if(window_size >= sequence_full_length) {
    starts <- 1
    ends <- sequence_full_length
  } else {
    starts <- genomic_bins_starts(start = 1, end = sequence_full_length, bin_size = window_size)
    starts <- starts[starts < sequence_full_length]
    if (length(starts) == 1) {
      ends <- sequence_full_length
    } else {
      ends <- c((starts[2 : length(starts)] - 1), sequence_full_length) + window_size # This makes overlapping windows!
    }
    ends[ends > sequence_full_length] <- sequence_full_length
  }

  # Divide into chunks
  scores <- NULL
  wins_per_chunk <- 100
  chunk_starts <- seq(1, length(starts), wins_per_chunk)
  chunk_starts <- c(chunk_starts, (length(starts) + 1))
  cat("Sequence full length: ", round(sequence_full_length / 1000000, 3), " Mbp \t", sep = "")
  cat("Chunks to complete: ", (length(chunk_starts) - 1), ". Finished: ", sep = "")
  date <- Sys.Date()
  for(i in 1 : (length(chunk_starts) - 1)) {
    sequence_substring <- fasta_sequence[starts[chunk_starts[i]] : ends[chunk_starts[i + 1] - 1]]
    for (j in (chunk_starts[i] : (chunk_starts[i + 1] - 1))) {
      result <- seq_win_score_int(1, window_size, kmer, sequence_substring[(starts[j] - starts[chunk_starts[i]] + 1) : (ends[j] - starts[chunk_starts[i]] + 1)], fraction_p)
      save(result, file = paste0(output_dir, "/", i, "_", j, "_", date, "_sequence_window_score_data"))
      remove(result)
    }
    cat(i, "")
    for(j in (chunk_starts[i] : (chunk_starts[i + 1] - 1))) {
      load(paste0(output_dir, "/", i, "_", j, "_", date, "_sequence_window_score_data"))
      unlink(paste0(output_dir, "/", i, "_", j, "_", date, "_sequence_window_score_data"))
      scores <- c(scores, result)
      remove(result)
    }
  }

  cat("\n")

  if ((sum(is.na(scores))>0) || (length(scores) == 0)) {
    stop("sequence_window_score did not produce valid result")
  }
  return(scores)
}
