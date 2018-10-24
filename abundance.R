seq_abundance <- function(seqs){
  abundance <- new.env()
  for (s in seqs){
    sequence <- s$sequence
    if (exists(sequence, envir = abundance)){
      abundance[[sequence]] <-abundance[[sequence]] + 1 
    }
    else{
      abundance[[sequence]] <- 0
    }
  }
  
  abundance
}
