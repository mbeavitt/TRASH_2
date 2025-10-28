main <- function(cmd_arguments) {
  # Helper functions
  log_hash <- function() {
    cat(strrep("#", 80), "\n", sep = "")
  }

  log_sep <- function() {
    cat(strrep("=", 80), "\n", sep = "")
  }

  set.seed(42)

  log_step <- function(step, total, description) {
    cat(sprintf("\n### %02d / %02d %s ### %s\n", step, total, description, Sys.time()))
    log_hash()
  }

  track_time <- function(event, data_type = "none", data_value = 0) {
    times$time <<- append(times$time, as.numeric(Sys.time()))
    times$event <<- append(times$event, event)
    times$data_type <<- append(times$data_type, data_type)
    times$data_value <<- append(times$data_value, data_value)
  }

  make_output_path <- function(suffix) {
    file.path(cmd_arguments$output_folder,
              paste0(basename(cmd_arguments$fasta_file), suffix))
  }

  make_temp_path <- function(...) {
    file.path(cmd_arguments$output_folder, paste(..., sep = "_"))
  }

  # Initialize
  log_hash()
  cat("### TRASH: workspace initialised ###", Sys.time(), "\n")
  log_hash()
  log_sep()

  mafft_dir <- "mafft"
  nhmmer_dir <- "nhmmer"
  log_messages <- ""

  # Settings
  kmer <- 10
  window_size <- round((cmd_arguments$max_rep_size + kmer) * 1.1)
  report_runtime <- TRUE
  add_sequence_info <- TRUE

  times <- list(time = as.numeric(Sys.time()), event = "Start main function", data_type = "none", data_value = 0)

  # Load fasta
  log_step(3, 13, paste("Loading the fasta file:", basename(cmd_arguments$fasta_file)))
  fasta_content <- load_and_validate_fasta(cmd_arguments$fasta_file)
  if (is.null(fasta_content)) {
    return(1)
  }
  log_sep()
  track_time("03 Fasta loaded", "Fasta total length", sum(sapply(fasta_content, length)))

  # Score sequences and identify repetitive regions
  log_step(4, 13, "Scoring sequences for repeats and identifying regions")
  repetitive_regions <- score_sequences_for_repeats(fasta_content, window_size, kmer,
                                                    cmd_arguments$output_folder, log_messages)
  if (is.null(repetitive_regions)) {
    print("No regions with repeats identified")
    return(0)
  }
  write.csv(repetitive_regions, make_output_path("_regarrays.csv"), row.names = FALSE)
  log_sep()
  track_time("Finished scoring and region identification", "Total region length",
             sum(repetitive_regions$ends - repetitive_regions$starts))
  # 06 / 14 Split regions into arrays 
  log_step(6, 13, "Identifying individual arrays with repeats")
  date <- Sys.Date()
  regions_per_chunk <- 100
  arrays <- NULL
  for (i in seq_along(fasta_content)) {
    cat("  Fasta sequence ", i, ": ", names(fasta_content)[i], " \t", sep = "")
    repetitive_regions_chr <- repetitive_regions[repetitive_regions$numID == i,]
    cat("Repetitive regions in the sequence: ", length(repetitive_regions_chr), " \t", sep = "")
    if(nrow(repetitive_regions_chr) == 0) {
      cat("\n")
      next
    }
    region_chunk <- seq(1, nrow(repetitive_regions_chr), regions_per_chunk) # divide into up to 100 data frame entries chunks on each chromosome, so up to 100 parallel, too much of a fasta is not good sent into the parallel
    cat("Chunks to complete: ", length(region_chunk), ". Finished: ", sep = "")
    region_chunk <- c(region_chunk, (nrow(repetitive_regions_chr) + 1))
    for(j in 1 : (length(region_chunk) - 1)) {
      sequence_substring <- fasta_content[[i]][repetitive_regions_chr$starts[region_chunk[j]] : (repetitive_regions_chr$ends[region_chunk[j + 1] - 1])]
      start_adjust <- repetitive_regions_chr$starts[region_chunk[j]] - 1
      for (k in (region_chunk[j] : (region_chunk[j+1] - 1))) {
        out <- split_and_check_arrays(start = repetitive_regions_chr$starts[k],
                                      end = repetitive_regions_chr$ends[k],
                                      sequence = sequence_substring[(repetitive_regions_chr$starts[k] - start_adjust) : (repetitive_regions_chr$ends[k] - start_adjust)],
                                      seqID = repetitive_regions_chr$seqID[k],
                                      numID = repetitive_regions_chr$numID[k],
                                      arrID = k,
                                      max_repeat = cmd_arguments$max_rep_size,
                                      min_repeat = cmd_arguments$min_rep_size,
                                      mafft = mafft_dir,
                                      temp_dir = cmd_arguments$output_folder,
                                      src_dir = getwd(),
                                      sink_output = FALSE,
                                      kmer = kmer)
        save(out, file = make_temp_path(i, j, k, date, "06_data"))
        remove(out)
      }
      cat(j, "")
      for (k in (region_chunk[j] : (region_chunk[j+1] - 1))) {
        temp_file <- make_temp_path(i, j, k, date, "06_data")
        load(temp_file)
        unlink(temp_file)
        arrays <- rbind(arrays, out)
        remove(out)
      }
    }
    cat("\n")
  }
  remove(repetitive_regions)
  log_sep()
  write.csv(arrays, make_output_path("_aregarrays.csv"), row.names = FALSE)

  track_time("Finished 06 split regions into arrays", "Total length of arrays",
             sum(arrays$end - arrays$start))
  # 07 / 14 Shift representative repeats and apply templates 
  log_step(7, 13, "Shifting representative and comparing templates")
  pb <- txtProgressBar(min = 0, max = nrow(arrays), style = 1)
  if (cmd_arguments$templates != 0) {
    templates <- read_fasta_and_list(cmd_arguments$templates)
    length_templates <- length(templates)
    if (length(templates) == 0) stop("No templates found within the template file")
  } else {
    templates <- 0
    length_templates <- 0
  }
  track_time("Finished 07 get templates", "Number of templates", length_templates)

  shifted_representatives <- character(nrow(arrays))
  for (i in seq_len(nrow(arrays))) {
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
    if (!inherits(arrays$representative[i], "character")) {
      shifted_representatives[i] <- "_"
    } else {
      shifted_representatives[i] <- shift_and_compare(arrays$representative[i], templates)
    }
  }
  arrays$representative <- shifted_representatives
  track_time("Finished 07 shift representatives and apply templates",
             "Number of templates times nrow arrays", length_templates * nrow(arrays))

  arrays$class <- ""
  for (i in seq_len(nrow(arrays))) {
    if (arrays$representative[i] == "_split_") {
      arrays$representative[i] <- ""
      next()
    }
    arrays$class[i] <- strsplit(arrays$representative[i], split = "_split_")[[1]][1]
    arrays$representative[i] <- strsplit(arrays$representative[i], split = "_split_")[[1]][2]
  }
  close(pb)
  track_time("Finished 07", "Arrays nrow", nrow(arrays))

  # 08 / 14 Classify unclassified and shift 
  log_step(8, 13, "Classifying remaining representative repeats")
  date <- Sys.Date()
  arrays <- classify_repeats(repeat_df = arrays)
  track_time("Finished 08 classify repeats", "Arrays nrow", nrow(arrays))

  write.csv(arrays, make_output_path("_arrays.csv"), row.names = FALSE)

  classes <- unique(arrays$class)
  classes <- classes[!(classes %in% c(names(templates), "none_identified"))]
  if (length(classes) != 0) {
    pb <- txtProgressBar(style = 1, min = 0, max = length(classes))
    for (i in seq_along(classes)) {
      arrays_class <- arrays[arrays$class == classes[i], ]
      arrays_class$representative <- shift_classes(arrays_class, kmer = 6)
      setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      temp_file <- make_temp_path(i, date, "08_data")
      save(arrays_class, file = temp_file)
      remove(arrays_class)
    }
    arrays_t <- NULL
    for (i in seq_along(classes)) {
      temp_file <- make_temp_path(i, date, "08_data")
      load(temp_file)
      unlink(temp_file)
      arrays_t <- rbind(arrays_t, arrays_class)
      remove(arrays_class)
    }
    arrays <- rbind(arrays_t,
                    arrays[which(arrays$class %in% c(names(templates), "none_identified")), ])
    close(pb)
    remove(arrays_t)
  } else {
    log_sep()
  }
  arrays <- arrays[order(arrays$start), ]
  arrays <- arrays[order(arrays$seqID), ]
  arrays$array_num_ID <- seq_len(nrow(arrays))

  arrays_no_representative <- arrays[arrays$class == "none_identified", ]
  arrays <- arrays[arrays$class != "none_identified", ]
  write.csv(arrays_no_representative, make_output_path("_no_repeats_arrays.csv"),
            row.names = FALSE)
  write.csv(arrays, make_output_path("_classarrays.csv"), row.names = FALSE)

  track_time("Finished 08 shift classes", "Unique classes number",
             length(unique(arrays$class)))

  remove(arrays_no_representative)

  if(nrow(arrays) == 0) {
    cat("No arrays with tandem repeats found under the settings\n")
    return(0)
  }

  # 09 / 14 Map repeats 
  log_step(9, 13, "Mapping array representatives")

  repeats <- NULL
  default_df = data.frame(seqID = vector(mode = "character"),
                          arrayID = vector(mode = "numeric"),
                          start = vector(mode = "numeric"),
                          end = vector(mode = "numeric"),
                          strand = vector(mode = "character"),
                          score = vector(mode = "numeric"),
                          eval = vector(mode = "numeric"),
                          width = vector(mode = "numeric"),
                          class = vector(mode = "character"),
                          representative = vector(mode = "character"),
                          score_template = vector(mode = "numeric"))
  arrays_per_chunk <- 100

  for(chromosome in seq_along(fasta_content)) {
    # For each chromosome 
    cat("  Fasta sequence ", chromosome, ": ", names(fasta_content)[chromosome], " \t", sep = "")
    arrays_chr <- arrays[arrays$numID == chromosome,]
    cat("Arrays in the sequence: ", nrow(arrays_chr), " \t", sep = "")
    if(nrow(arrays_chr) == 0) {
      cat("\n")
      next
    }
    array_chunk <- seq(1, nrow(arrays_chr), arrays_per_chunk) # divide into up to 100 data frame entries chunks on each chromosome, so up to 100 parallel, too much of a fasta is not good sent into the parallel
    cat("Chunks to complete: ", length(array_chunk), ". Finished: ", sep = "")
    array_chunk <- c(array_chunk, (nrow(arrays_chr) + 1))
    for(j in 1 : (length(array_chunk) - 1)) {
      sequence_substring <- fasta_content[[chromosome]][arrays_chr$start[array_chunk[j]] : (arrays_chr$end[array_chunk[j + 1] - 1])]
      adjust_start <- arrays_chr$start[array_chunk[j]] - 1
      arrays_chunk_IDs <- array_chunk[j] : (array_chunk[j + 1] - 1)
      for (i in arrays_chunk_IDs) {
        if (arrays_chr$representative[i] == "") {
          cat(i, "")
          next
        }
        array_sequence <- sequence_substring[(arrays_chr$start[i] - adjust_start) : (arrays_chr$end[i] - adjust_start)]
        cat(i, "_ ", sep = "")
        if (arrays_chr$top_N[i] >= 14) {
          # nhmmer for repeats of 14+ bp
          repeats_df <- map_nhmmer(cmd_arguments$output_folder,
                                   arrayID = arrays_chr$array_num_ID[i],
                                   arrays_chr$representative[i],
                                   arrays_chr$seqID[i],
                                   arrays_chr$start[i],
                                   arrays_chr$end[i],
                                   array_sequence,
                                   nhmmer_dir)
        } else {
          # matchpattern for shorter
          repeats_df <- map_default(arrayID = arrays_chr$array_num_ID[i],
                                    arrays_chr$representative[i],
                                    arrays_chr$seqID[i],
                                    arrays_chr$start[i],
                                    paste(array_sequence, collapse = ""))
        }
        if (nrow(repeats_df) < 2) {
          cat(i, "")
          next
        }
        # add width and class 
        repeats_df$width <- repeats_df$end - repeats_df$start + 1
        repeats_df$class <- arrays_chr$class[i]
        # Handle overlaps
        repeats_df <- handle_overlaps(repeats_df, overlap_threshold = 0.1)
        if (nrow(repeats_df) < 3) {
          cat(i, "")
          next
        }
        # handle gaps if proper array
        repeats_df <- handle_gaps(repeats_df, representative_len = arrays_chr$top_N[i])
        # Skip if handle_gaps removed all repeats
        if (nrow(repeats_df) < 3) {
          cat(i, "")
          next
        }
        # Recalculate representative 
        # TODO: Use more repeats for long arrays, this is stringent and good for most arrays, but some deserve a better recalculation
        max_repeats_to_align <- 15
        min_repeats_to_recalculate <- 10
        repeats_df$representative <- arrays_chr$representative[i]
        repeats_df$score_template <- -1
        sample_IDs <- which(repeats_df$strand != ".")
        if (length(sample_IDs) >= min_repeats_to_recalculate) {
          if (length(sample_IDs) > max_repeats_to_align) {
            sample_IDs <- sample(sample_IDs, max_repeats_to_align)
          }
          repeats_seq <- unlist(lapply(sample_IDs, function(X) {
            paste0(sequence_substring[(repeats_df$start[X] - adjust_start):
                                      (repeats_df$end[X] - adjust_start)],
                   collapse = "")
          }))
          strands <- repeats_df$strand[sample_IDs]
          repeats_seq[which(strands == "-")] <- unlist(
            lapply(repeats_seq[which(strands == "-")], rev_comp_string))
          alignment <- write_align_read(
            mafft_exe = mafft_dir,
            temp_dir = cmd_arguments$output_folder,
            sequences = repeats_seq,
            name = paste(basename(cmd_arguments$fasta_file),
                        arrays_chr$seqID[i],
                        arrays_chr$array_num_ID[i],
                        i, runif(1, 0, 1), sep = "_"))
          consensus <- consensus_N(alignment, arrays_chr$top_N[i])
          if (length(consensus) != 0) repeats_df$representative <- consensus
          remove(alignment, consensus, repeats_seq, strands)
        }
        # check if short gaps contain the repeat 
        repeats_df <- fill_gaps(repeats_df, array_sequence, arrays_chr$start[i])
        # Change to edit distance based score
        repeats_seq <- unlist(lapply(seq_len(nrow(repeats_df)), function(X) {
          paste0(sequence_substring[(repeats_df$start[X] - adjust_start):
                                    (repeats_df$end[X] - adjust_start)],
                 collapse = "")
        }))
        costs <- list(insertions = 1, deletions = 1, substitutions = 1)
        rep_len <- nchar(repeats_df$representative[1])
        plus_strand <- repeats_df$strand == "+"
        minus_strand <- repeats_df$strand == "-"

        if (sum(plus_strand) > 0) {
          repeats_df$score[plus_strand] <- adist(repeats_df$representative[1],
                                                  repeats_seq[plus_strand],
                                                  costs)[1, ] / rep_len * 100
        }
        if (sum(minus_strand) > 0) {
          repeats_df$score[minus_strand] <- adist(rev_comp_string(repeats_df$representative[1]),
                                                   repeats_seq[minus_strand])[1, ] / rep_len * 100
        }
        if (arrays_chr$class[i] %in% names(templates)) {
          template <- paste(templates[[which(names(templates) == arrays_chr$class[i])]], collapse = "")
          temp_len <- nchar(template)
          if (sum(plus_strand) > 0) {
            repeats_df$score_template[plus_strand] <- adist(template,
                                                            repeats_seq[plus_strand],
                                                            costs)[1, ] / temp_len * 100
          }
          if (sum(minus_strand) > 0) {
            repeats_df$score_template[minus_strand] <- adist(rev_comp_string(template),
                                                             repeats_seq[minus_strand])[1, ] / temp_len * 100
          }
        }
        # Correct repeats split into two by nhmmer
        score_min_to_merge <- 30
        size_max_to_merge <- 1.0
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
              new_score <- adist(rep_to_compare, merged_seq, costs)[1, ] / rep_len * 100

              if (new_score < min(repeats_df$score[i_r:(i_r + 1)])) {
                repeats_df$end[i_r] <- repeats_df$end[i_r + 1]
                repeats_df <- repeats_df[-(i_r + 1), ]
                repeats_seq <- repeats_seq[-(i_r + 1)]
                repeats_seq[i_r] <- paste0(sequence_substring[(repeats_df$start[i_r] - adjust_start):
                                                              (repeats_df$end[i_r] - adjust_start)],
                                          collapse = "")
                repeats_df$score[i_r] <- new_score
                if (arrays_chr$class[i_r] %in% names(templates)) {
                  template <- paste(templates[[which(names(templates) == arrays_chr$class[i_r])]], collapse = "")
                  temp_to_compare <- if (both_plus) template else rev_comp_string(template)
                  repeats_df$score_template[i_r] <- adist(temp_to_compare, repeats_seq[i_r])[1, ] /
                                                    nchar(template) * 100
                }
              }
            }
          }
          i_r <- i_r + 1
        }
        # Handle edge repeats
        repeats_df <- handle_edge_repeat(repeats_df, sequence_substring, adjust_start)
        # TODO: make sure the edge repeats have their template score recalculated too

        remove(repeats_seq)
        if (nrow(repeats_df) < 3) {
          cat(i, "")
          next
        }
        repeats_df <- repeats_df[c("seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class", "representative", "score_template")]
        repeats_df$seqID <- as.character(repeats_df$seqID)
        repeats_df$arrayID <- as.numeric(repeats_df$arrayID)
        repeats_df$start <- as.numeric(repeats_df$start)
        repeats_df$end <- as.numeric(repeats_df$end)
        repeats_df$strand <- as.character(repeats_df$strand)
        repeats_df$score <- as.numeric(repeats_df$score)
        repeats_df$eval <- as.numeric(repeats_df$eval)
        repeats_df$width <- as.numeric(repeats_df$width)
        repeats_df$class <- as.character(repeats_df$class)
        repeats_df$representative <- as.character(repeats_df$representative)
        repeats_df$score_template <- as.numeric(repeats_df$score_template)
        save(repeats_df, file = make_temp_path(i, date, "09_data"))
        remove(repeats_df)
        cat(i, "")
      }
      for(i in arrays_chunk_IDs) {
        temp_file <- make_temp_path(i, date, "09_data")
        if(file.exists(temp_file)) {
          load(temp_file)
          unlink(temp_file)
          if(sum(names(repeats) != names(repeats_df)) > 0) {
          print(paste(i, names(repeats), names(repeats_df)))
          print(str(repeats))
          print(str(repeats_df))
        }
        if(nrow(repeats_df) > 0) repeats <- rbind(repeats, repeats_df)
        remove(repeats_df)
        }
      }
    }
    cat("\n")
  }

  track_time("Finished 09 map repeats", "Repeats number", nrow(repeats))

  # 11 / 14 Summarise array information (where's 10 / 14??)
  log_step(11, 13, "Summarising array information")

  for (i in seq_len(nrow(arrays))) {
    if (sum((repeats$arrayID == arrays$array_num_ID[i]) > 0)) {
      arrays$representative[i] <- repeats$representative[which(repeats$arrayID == i)[1]]
    }
  }
  repeats <- repeats[c("seqID", "arrayID", "start", "end", "strand", "score",
                       "eval", "width", "class", "score_template")]

  track_time("Finished 11 reassign array representatives", "Arrays nrow", nrow(arrays))

  arrays$repeats_number <- 0
  arrays$median_repeat_width <- 0
  arrays$median_score <- -1

  for (i in seq_len(nrow(arrays))) {
    repeats_temp <- repeats[repeats$arrayID == arrays$array_num_ID[i], ]
    if (nrow(repeats_temp) > 0) {
      arrays$repeats_number[i] <- nrow(repeats_temp)
      arrays$median_repeat_width[i] <- ceiling(median(repeats_temp$width))
      arrays$median_score[i] <- ceiling(median(repeats_temp$score))
    }
  }
  arrays <- arrays[arrays$repeats_number != 0, ]
  log_sep()

  track_time("Finished 11 summarise array info", "Arrays nrow", nrow(arrays))

  # 12 / 14 Save array output
  log_step(12, 13, "Saving the array table")
  write.csv(arrays, make_output_path("_arrays.csv"), row.names = FALSE)
  export_gff(annotations.data.frame = arrays,
             output = cmd_arguments$output_folder,
             file.name = paste0(basename(cmd_arguments$fasta_file), "_arrays"),
             source = "TRASH",
             type = "Satellite_array",
             seqid = 3,
             start = 1,
             end = 2,
             score = 5,
             attributes = c(9, 10, 11),
             attribute.names = c("Name=", "Repeat_no=", "Repeat_median_width="))
  log_sep()
  track_time("Finished 12 saved array info", "Arrays nrow", nrow(arrays))

  # 13 / 14 Save repeat output
  log_step(13, 13, "Saving the repeats table")
  write.csv(repeats, make_output_path("_repeats.csv"), row.names = FALSE)
  export_gff(annotations.data.frame = repeats,
             output = cmd_arguments$output_folder,
             file.name = paste0(basename(cmd_arguments$fasta_file), "_repeats"),
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

  log_sep()
  track_time("Finished 13 saved repeats info", "Repeats number", nrow(repeats))

  # 14 / 14 Done
  if (report_runtime) {
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
}
