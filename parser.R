convert_quality <- function(line) {
  # TODO: accommodate other quality score conversions (e.g., Solexa)
  # @arg line:
  # @return <describe return value>
  
  # remove trailing whitespace (incl. line break)
  line <- strsplit(sub("\\s+$", "", line), "")[[1]]
  
  result <- vector()
  for (letter in line){
    score <- utf8ToInt(letter) - 33
    if (score < 0 | score > 41){
      stop(paste("ERROR: Unexpected integer value in convert_quality():", score))
    } else{
      result <- c(result, score)
    }
  }
  result
}

parse_fastq <- function(file) {
  # @arg file: 
  con <- file(file, "r")
  # con <- file
  tmp.list <- list()
  header <- NULL
  sequence <- ""
  quality <- vector()
  ln <- 0
  seqs <- 0
  
  while( TRUE ) {
    line = readLines(con, n=1)
    if (length(line) == 0) {
      break  # reached end of file
    }
    position <- ln %% 4
    
    if (position == 0 && startsWith(line, "@")) {
      if (!is.null(header)){
        seqs <- seqs + 1
        tmp.list[[seqs]] <- list(header, sequence, quality)
      }
      header <- sub("^@", "", line)
    }
    else if (position == 1) {
      sequence <- line
    }
    else if (position == 2 && startsWith(line, "+")) {
      ln <- ln + 1
      next
    }
    else if (position == 3){
      quality <- convert_quality(line)
    }
    else{
      stop(paste("ERROR: Failed to parse FASTQ at line:\n", line))
    }
    ln <- ln + 1
  }
  close(con)
  
  # handle last entry
  seqs <- seqs + 1
  tmp.list[[seqs]] <- list(header, sequence, quality)
  tmp.list
}

