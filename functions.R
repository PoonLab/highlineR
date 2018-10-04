import.files <- function(path, pat, ...){
  tmp.list.1 <- list.files(path, pattern = pat)
  tmp.list.2 <- list(length = length(tmp.list.1))
  for (i in 1:length(tmp.list.1)){
    tmp.list.2[[i]] <- readLines(paste0(path, tmp.list.1[i]))
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

seq_abundance <- function(seqs, ...){
  no_desc <- seqs[c(FALSE, TRUE)]
  abundance <- new.env()
  for (s in no_desc){
    if (exists(s, envir = abundance)){
      abundance[[s]] <-abundance[[s]] + 1 
    }
    else{
      abundance[[s]] <- 0
    }
  }
  abundance
}

# trial
fastq.import <- import.fastq("/home/lisamonique/Documents/highlineR/fastq/")
fasta <- fastqtofasta(fastq.import[[1]])
abundance <- seq_abundance(fasta)
ls.str(abundance)
