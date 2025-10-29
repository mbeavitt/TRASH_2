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

  mafft_executable <- "mafft"
  nhmmer_executable <- "nhmmer"
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
                                      mafft = mafft_executable,
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
                                   nhmmer_executable)
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
        array_info <- list(
          representative = arrays_chr$representative[i],
          seqID = arrays_chr$seqID[i],
          array_num_ID = arrays_chr$array_num_ID[i],
          top_N = arrays_chr$top_N[i],
          i = i
        )
        repeats_df <- recalculate_representative(
          repeats_df, sequence_substring, adjust_start, array_info,
          mafft_executable, cmd_arguments$output_folder, basename(cmd_arguments$fasta_file)
        )
        # Check if short gaps contain the repeat
        repeats_df <- fill_gaps(repeats_df, array_sequence, arrays_chr$start[i])

        # Calculate edit distance scores
        score_result <- calculate_repeat_edit_distance_scores(
          repeats_df, sequence_substring, adjust_start, templates, arrays_chr$class[i]
        )
        repeats_df <- score_result$repeats_df
        repeats_seq <- score_result$repeats_seq
        # Merge repeats that were incorrectly split by nhmmer
        merge_result <- merge_split_repeats(
          repeats_df, repeats_seq, sequence_substring, adjust_start,
          templates, arrays_chr$class[i]
        )
        repeats_df <- merge_result$repeats_df
        repeats_seq <- merge_result$repeats_seq
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

  # Save outputs
  log_step(12, 13, "Saving outputs")
  save_array_output(arrays, cmd_arguments$output_folder,
                   basename(cmd_arguments$fasta_file), make_output_path)
  log_sep()
  track_time("Saved array output", "Arrays nrow", nrow(arrays))

  save_repeat_output(repeats, fasta_content, cmd_arguments$output_folder,
                    basename(cmd_arguments$fasta_file), make_output_path, add_sequence_info)
  log_sep()
  track_time("Saved repeat output", "Repeats number", nrow(repeats))

  # Generate runtime report
  if (report_runtime) {
    generate_runtime_report(times, make_output_path)
  }
}
