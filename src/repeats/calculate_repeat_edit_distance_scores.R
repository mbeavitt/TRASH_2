calculate_repeat_edit_distance_scores <- function(repeats_df, sequence_substring,
                                                   adjust_start, templates, array_class) {
  # Extract sequences for all repeats
  repeats_seq <- unlist(lapply(seq_len(nrow(repeats_df)), function(X) {
    paste0(sequence_substring[(repeats_df$start[X] - adjust_start):
                              (repeats_df$end[X] - adjust_start)],
           collapse = "")
  }))

  # Calculate edit distance scores
  plus_strand <- repeats_df$strand == "+"
  minus_strand <- repeats_df$strand == "-"

  # Score against representative
  if (sum(plus_strand) > 0) {
    repeats_df$score[plus_strand] <- calculate_edit_distance_score(
      repeats_df$representative[1],
      repeats_seq[plus_strand]
    )
  }
  if (sum(minus_strand) > 0) {
    repeats_df$score[minus_strand] <- calculate_edit_distance_score(
      rev_comp_string(repeats_df$representative[1]),
      repeats_seq[minus_strand]
    )
  }

  # Score against template if applicable
  if (array_class %in% names(templates)) {
    template <- paste(templates[[which(names(templates) == array_class)]], collapse = "")
    if (sum(plus_strand) > 0) {
      repeats_df$score_template[plus_strand] <- calculate_edit_distance_score(
        template,
        repeats_seq[plus_strand]
      )
    }
    if (sum(minus_strand) > 0) {
      repeats_df$score_template[minus_strand] <- calculate_edit_distance_score(
        rev_comp_string(template),
        repeats_seq[minus_strand]
      )
    }
  }

  return(list(repeats_df = repeats_df, repeats_seq = repeats_seq))
}
