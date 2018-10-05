import.files <- function(path, pat, ...){
  tmp.list.1 <- list.files(path, pattern = pat)
  tmp.list.2 <- list(length = length(tmp.list.1))
  for (i in 1:length(tmp.list.1)){
    # tmp.list.2[[i]] <- readLines(paste0(path, tmp.list.1[i]))
    tmp.list.2[[i]] <- file(paste0(path, tmp.list.1[i]), "r")
  }
  names(tmp.list.2) <- tmp.list.1
  tmp.list.2
}

import.fastq <- function(path, ...){
  import.files(path, ".(fastq|fq)$")
}

fastqtofasta <- function(seqs, ...){
  seqs[c(TRUE, TRUE, FALSE, FALSE)]
}

