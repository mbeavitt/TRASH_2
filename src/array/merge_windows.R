merge_windows <- function(list_of_scores, window_size, sequence_full_length) {
  # TODO make the treshold dynamic
  threshold <- 90

  # if(sequence_full_length < window_size) {
  #   return(data.frame(starts = vector(mode = "numeric"),
  #                     ends = vector(mode = "numeric"),
  #                     scores = vector(mode = "numeric")))
  # }

  if (sum(list_of_scores < threshold) == 0) {
    return(data.frame(
                      starts = vector(mode = "numeric"),
                      ends = vector(mode = "numeric"),
                      scores = vector(mode = "numeric")
    ))
  }
  result <- genomic_bins(start = 1, end = sequence_full_length, bin_size = window_size)
  starts <- result$starts
  ends <- result$ends

  if (length(starts) == 1) {
    repetitive_regions <- data.frame(
                                     starts = starts,
                                     ends = sequence_full_length,
                                     scores = list_of_scores
    )
    return(repetitive_regions)
  }

  if(length(starts) != length(ends)) {
    stop("merge_windows starts != ends")
  }
  if(length(starts) != length(list_of_scores)) {
    stop("merge_windows starts != list_of_scores")
  }
  if(length(ends) != length(list_of_scores)) {
    stop("merge_windows ends != list_of_scores")
  }
  
  repetitive_regions <- data.frame(
                                   starts = starts[list_of_scores < threshold],
                                   ends = ends[list_of_scores < threshold],
                                   scores = list_of_scores[list_of_scores < threshold]
  )
  
  if (nrow(repetitive_regions) < 2) return(repetitive_regions)

  i <- 1
  while (i < nrow(repetitive_regions)) {
    if ((repetitive_regions$ends[i] + 1) >= (repetitive_regions$starts[i + 1])) {
      repetitive_regions$scores[i] <- (
                                         (repetitive_regions$scores[i] * (repetitive_regions$ends[i] - repetitive_regions$starts[i])) +
                                         (repetitive_regions$scores[i + 1] * (repetitive_regions$ends[i + 1] - repetitive_regions$starts[i + 1]))
      ) / (
                                         repetitive_regions$ends[i + 1] - repetitive_regions$starts[i + 1] + repetitive_regions$ends[i] - repetitive_regions$starts[i]
      )
      repetitive_regions$ends[i] <- repetitive_regions$ends[i + 1]
      repetitive_regions <- repetitive_regions[-(i + 1), ]
      i <- i - 1
    }
    i <- i + 1
  }
  return(repetitive_regions)
}
