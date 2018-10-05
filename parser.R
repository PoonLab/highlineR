parse_fastq <- function(file){
  tmp.list <- list()
  header <- NULL
  sequence <- ""
  quality <- vector()
  ln = 0
  for (line in file){
    position <- ln %% 4
    
    if (position == 0 && startsWith(line, "@")){
      if (! is.null(header)){
        tmp.list <- c(tmp.list, c(header, sequence, quality))
      }
      else{
        header <- sub("\\s+$", "", line)
      }
    }
    else if (position == 1){
      sequence = sub("\\s+$", "", line)
    }
    else if (position == 2 && startsWith(line, "+")){
      ln = ln + 1
      next
    }
    else if (position == 3){
      quality <- convert_quality(line)
    }
    else{
      stop(paste("ERROR: Failed to parse FASTQ at line:\n", line))
    }
    ln = ln + 1
  }
  tmp.list <- c(tmp.list, c(header, sequence, quality))
  tmp.list
}