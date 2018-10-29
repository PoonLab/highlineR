seq_abundance <- function(seqs){
  # @arg seqs: list of sequences
  # @return environment of (sequence: abunance) key:value pairs
  
  abundance <- new.env()
  for (s in seqs){
    sequence <- s$sequence
    if (exists(sequence, envir = abundance)){
      # if sequence already in structure, increment count
      abundance[[sequence]] <-abundance[[sequence]] + 1 
    }
    else{
      # otherwise, add sequence and initiate count
      abundance[[sequence]] <- 0
    }
  }
  
  abundance
}
