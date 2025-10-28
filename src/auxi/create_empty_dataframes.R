create_empty_repeat_df <- function() {
  # Create an empty dataframe with standard repeat columns
  #
  # Returns:
  #   Empty dataframe with repeat structure

  data.frame(
    seqID = character(),
    arrayID = numeric(),
    start = numeric(),
    end = numeric(),
    strand = character(),
    score = numeric(),
    eval = numeric(),
    width = numeric(),
    class = character(),
    representative = character(),
    score_template = numeric(),
    stringsAsFactors = FALSE
  )
}

create_empty_array_df <- function() {
  # Create an empty dataframe with standard array columns
  #
  # Returns:
  #   Empty dataframe with array structure

  data.frame(
    start = numeric(),
    end = numeric(),
    top_N = numeric(),
    representative = character(),
    seqID = character(),
    numID = numeric(),
    stringsAsFactors = FALSE
  )
}
