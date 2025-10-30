seq_win_score_int <- function(start, end, kmer, fasta_extraction, fraction_p) {
  if((end - start) <= kmer) return(100)
  # kmers <- unlist(lapply(X = start : (end - kmer), FUN = extract_kmers, kmer, fasta_extraction))
  kmers <- unlist(lapply(X = (start : (end - kmer)), function(X) return(paste(fasta_extraction[X : (X + kmer - 1)], collapse = ""))))
  counts_kmers <- table(kmers)
  counts_kmers <- counts_kmers[!grepl("n", names(counts_kmers))]
  counts_kmers <- counts_kmers[!grepl("N", names(counts_kmers))]

  if (sum(counts_kmers) < (kmer * 2)) return(100)
  score <- (100 * sum(counts_kmers[counts_kmers == 1]) / sum(counts_kmers))
  remove(start, end, kmer, fasta_extraction, fraction_p, counts_kmers, kmers)
  gc()
  return(score)
}

