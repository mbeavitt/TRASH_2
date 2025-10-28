extract_sequence_string <- function(sequence_vector, start_pos, end_pos) {
  # Extract a substring from a sequence vector and collapse to a single string
  #
  # Args:
  #   sequence_vector: Vector of individual nucleotides
  #   start_pos: Starting position (1-indexed)
  #   end_pos: Ending position (inclusive)
  #
  # Returns:
  #   Single string containing the sequence

  paste0(sequence_vector[start_pos:end_pos], collapse = "")
}
