merge_split_repeats <- function(repeats_df, repeats_seq, sequence_substring,
                                adjust_start, templates, array_class,
                                score_min_to_merge = 30, size_max_to_merge = 1.0) {
  # Correct repeats split into two by nhmmer
  # Merge adjacent repeats if they have poor scores individually but good score combined

  i_r <- 1
  while(i_r < nrow(repeats_df)) {
    both_high_score <- repeats_df$score[i_r] > score_min_to_merge &&
                       repeats_df$score[i_r + 1] > score_min_to_merge
    combined_short <- sum(repeats_df$width[i_r:(i_r + 1)]) <
                     (size_max_to_merge * nchar(repeats_df$representative[1]))

    if (both_high_score && combined_short) {
      both_plus <- (repeats_df$strand[i_r] == "+") && (repeats_df$strand[i_r + 1] == "+")
      both_minus <- (repeats_df$strand[i_r] == "-") && (repeats_df$strand[i_r + 1] == "-")

      if (both_plus || both_minus) {
        merged_seq <- paste0(repeats_seq[i_r:(i_r + 1)], collapse = "")
        rep_to_compare <- if (both_plus) {
          repeats_df$representative[1]
        } else {
          rev_comp_string(repeats_df$representative[1])
        }
        new_score <- calculate_edit_distance_score(rep_to_compare, merged_seq)

        if (new_score < min(repeats_df$score[i_r:(i_r + 1)])) {
          # Merge the two repeats
          repeats_df$end[i_r] <- repeats_df$end[i_r + 1]
          repeats_df <- repeats_df[-(i_r + 1), ]
          repeats_seq <- repeats_seq[-(i_r + 1)]

          # Update sequence and score
          repeats_seq[i_r] <- extract_sequence_string(sequence_substring, (repeats_df$start[i_r] - adjust_start), (repeats_df$end[i_r] - adjust_start))
          repeats_df$score[i_r] <- new_score

          # Update template score if applicable
          if (array_class %in% names(templates)) {
            template <- paste(templates[[which(names(templates) == array_class)]], collapse = "")
            temp_to_compare <- if (both_plus) template else rev_comp_string(template)
            repeats_df$score_template[i_r] <- calculate_edit_distance_score(temp_to_compare, repeats_seq[i_r])
          }
        }
      }
    }
    i_r <- i_r + 1
  }

  return(list(repeats_df = repeats_df, repeats_seq = repeats_seq))
}
