recalculate_representative <- function(repeats_df, sequence_substring, adjust_start,
                                       array_info, mafft_dir, temp_dir, fasta_basename,
                                       max_repeats_to_align = 15, min_repeats_to_recalculate = 10) {
  # Initialize with existing representative
  repeats_df$representative <- array_info$representative
  repeats_df$score_template <- -1

  # Sample repeats for alignment
  sample_IDs <- which(repeats_df$strand != ".")
  if (length(sample_IDs) < min_repeats_to_recalculate) {
    return(repeats_df)
  }

  # Limit number of repeats to align
  if (length(sample_IDs) > max_repeats_to_align) {
    sample_IDs <- sample(sample_IDs, max_repeats_to_align)
  }

  # Extract sequences for selected repeats
  repeats_seq <- unlist(lapply(sample_IDs, function(X) {
    paste0(sequence_substring[(repeats_df$start[X] - adjust_start):
                              (repeats_df$end[X] - adjust_start)],
           collapse = "")
  }))

  # Reverse complement minus strand sequences
  strands <- repeats_df$strand[sample_IDs]
  repeats_seq[which(strands == "-")] <- unlist(
    lapply(repeats_seq[which(strands == "-")], rev_comp_string))

  # Perform multiple sequence alignment
  alignment <- write_align_read(
    mafft_exe = mafft_dir,
    temp_dir = temp_dir,
    sequences = repeats_seq,
    name = paste(fasta_basename, array_info$seqID, array_info$array_num_ID,
                array_info$i, runif(1, 0, 1), sep = "_"))

  # Generate consensus sequence
  consensus <- consensus_N(alignment, array_info$top_N)
  if (length(consensus) != 0) {
    repeats_df$representative <- consensus
  }

  remove(alignment, consensus, repeats_seq, strands)
  return(repeats_df)
}
