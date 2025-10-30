seq_win_score_int <- function(start, end, kmer, fasta_extraction) {
  if ((end - start) <= kmer) return(100)

  kmers <- sapply(start:(end - kmer), function(i)
    paste(fasta_extraction[i:(i + kmer - 1)], collapse = "")
  )

  counts <- table(kmers[!grepl("[nN]", kmers)])

  if (sum(counts) < (kmer * 2)) return(100)

  score <- 100 * sum(counts == 1) / sum(counts)
  return(score)
}

