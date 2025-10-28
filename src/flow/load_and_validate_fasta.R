load_and_validate_fasta <- function(fasta_file) {
  # Load fasta file
  fasta_content <- read_fasta_and_list(fasta_file)

  # Check if empty
  if (length(fasta_content) == 0) {
    warning("Fasta could not be read or is empty")
    return(NULL)
  }

  # Handle duplicate sequence names
  if (length(names(fasta_content)) != length(unique(names(fasta_content)))) {
    msg <- paste0("\nWARNING: Sequence names in the ", basename(fasta_file),
                  " fasta file are not unique \n They were appended to avoid assignment errors \n")
    warning(msg)
    cat(msg, "\n", "Adjustments made: \n", sep = "")

    fasta_names <- names(fasta_content)
    unique_names <- unique(fasta_names)
    for (i in seq_along(unique_names)) {
      matches <- fasta_names == unique_names[i]
      if (sum(matches) > 1) {
        old_names <- fasta_names[matches]
        new_names <- paste0(unique_names[i], seq_len(sum(matches)))
        cat("Old names:", old_names, "\n New names:", new_names, "\n\n\n")
        names(fasta_content)[matches] <- new_names
      }
    }
  }

  return(fasta_content)
}
