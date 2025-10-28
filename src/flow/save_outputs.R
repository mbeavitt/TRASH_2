save_array_output <- function(arrays, output_folder, fasta_basename, make_output_path) {
  cat("### Saving array output ###\n")
  write.csv(arrays, make_output_path("_arrays.csv"), row.names = FALSE)
  export_gff(annotations.data.frame = arrays,
             output = output_folder,
             file.name = paste0(fasta_basename, "_arrays"),
             source = "TRASH",
             type = "Satellite_array",
             seqid = 3,
             start = 1,
             end = 2,
             score = 5,
             attributes = c(9, 10, 11),
             attribute.names = c("Name=", "Repeat_no=", "Repeat_median_width="))
}

save_repeat_output <- function(repeats, fasta_content, output_folder, fasta_basename,
                                make_output_path, add_sequence_info = TRUE) {
  cat("### Saving repeat output ###\n")
  write.csv(repeats, make_output_path("_repeats.csv"), row.names = FALSE)
  export_gff(annotations.data.frame = repeats,
             output = output_folder,
             file.name = paste0(fasta_basename, "_repeats"),
             source = "TRASH",
             type = "Satellite_DNA",
             seqid = 1,
             start = 3,
             end = 4,
             strand = 5,
             attributes = c(9, 6, 10),
             attribute.names = c("Name=", "Arry_EDS=", "Family_EDS="))

  if (add_sequence_info) {
    repeats$sequence <- unlist(lapply(seq_len(nrow(repeats)), function(X) {
      paste0(fasta_content[[which(names(fasta_content) == repeats$seqID[X])]][
        repeats$start[X]:repeats$end[X]], collapse = "")
    }))

    minus_indices <- which(repeats$strand == "-")
    repeats$sequence[minus_indices] <- unlist(
      lapply(repeats$sequence[minus_indices], rev_comp_string))
    write.csv(repeats, make_output_path("_repeats_with_seq.csv"), row.names = FALSE)
  }
}

generate_runtime_report <- function(times, make_output_path) {
  times$time_passed <- 0
  times$time_per_Mbp <- 0
  times$time_per_event_data_value <- 0
  for (i in 2:length(times$time)) {
    times$time_passed <- append(times$time_passed, (times$time[i] - times$time[i - 1]))
    times$time_per_Mbp <- append(times$time_per_Mbp,
                                  (1000000 * times$time_passed[i] / times$data_value[2]))
    times$time_per_event_data_value <- append(times$time_per_event_data_value,
                                               (1000000 * times$time_passed[i] / times$data_value[i]))
  }
  write.csv(times, make_output_path("_run_time.csv"), row.names = FALSE)
}
