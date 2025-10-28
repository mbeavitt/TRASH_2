score_sequences_for_repeats <- function(fasta_content, window_size, kmer, output_folder, log_messages = "") {
  # Calculate repeat scores for each sequence
  cat("### Calculating repeat scores for each sequence ###\n")
  chromosome_lengths <- unlist(lapply(seq_along(fasta_content), function(X) length(fasta_content[[X]])))
  cat("  Assembly total length:\t", round(sum(chromosome_lengths) / 1000000, 1), "Mbp \n")
  cat("  Sequences count:\t\t\t", length(chromosome_lengths), " \n")
  cat("  Sequences names:\t\t\t", names(fasta_content), "\n")
  cat("  Sequences lengths (bp):\t", chromosome_lengths, "\n\n")

  repeat_scores <- list()
  for (i in seq_along(fasta_content)) {
    cat("  Fasta sequence ", i, ": ", names(fasta_content)[i], " \t", sep = "")
    repeat_scores <- append(repeat_scores,
                           list(sequence_window_score(fasta_content[[i]], window_size, kmer,
                                                      output_dir = output_folder)))
  }

  # Identify regions with high repeat content and merge into a dataframe
  cat("\n### Identifying regions with high repeat content ###\n")
  repetitive_regions <- data.frame(starts = NULL, ends = NULL, scores = NULL,
                                   seqID = NULL, numID = NULL)
  for (i in seq_along(repeat_scores)) {
    if (length(repeat_scores[[i]]) == 0) next
    regions_of_sequence <- merge_windows(list_of_scores = repeat_scores[[i]],
                                         window_size = window_size,
                                         sequence_full_length = length(fasta_content[[i]]),
                                         log_messages)
    if (nrow(regions_of_sequence) != 0) {
      regions_of_sequence$seqID <- names(fasta_content)[[i]]
      regions_of_sequence$numID <- i
      repetitive_regions <- rbind(repetitive_regions, regions_of_sequence)
    }
  }

  # Check if any regions were found
  if (!inherits(repetitive_regions, "data.frame") || nrow(repetitive_regions) == 0) {
    return(NULL)
  }

  return(repetitive_regions)
}
