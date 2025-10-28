calculate_edit_distance_score <- function(reference, target, costs = list(insertions = 1, deletions = 1, substitutions = 1)) {
  # Calculate edit distance score as a percentage
  #
  # Args:
  #   reference: Reference sequence string
  #   target: Target sequence(s) to compare (can be a vector)
  #   costs: List with insertions, deletions, and substitutions costs
  #
  # Returns:
  #   Edit distance as percentage of reference length (0-100)
  #   Lower scores indicate higher similarity

  ref_length <- nchar(reference)
  if (ref_length == 0) {
    return(100)  # Avoid division by zero
  }

  distances <- adist(reference, target, costs)[1, ]
  scores <- (distances / ref_length) * 100

  return(scores)
}
