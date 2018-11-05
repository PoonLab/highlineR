#TODO: verify parsed list empty before parsing

parse <- function(x) {
  UseMethod("parse", x)
}

parse.session <- function(session = highlineR.data) {
  #arg session: environment containing imported sequence Data objects
  
  for (data in ls(session)) {
    parse(get(data, envir = session, inherits = FALSE))
  }
}

parse.fasta <- function(data) {
  # @arg data: Data object containing absolute or relative path to a FASTA file
  # populates Data object's raw_seq list with (header, sequence) lists
  
  con <- file(data$path, "r")
  
  # prepare containers
  header <- NULL
  sequence <- ""
  
  while( TRUE ) {
    line = readLines(con, n = 1)
    if (length(line) == 0) {
      break  # reached end of file
    }
    
    if (startsWith(line, ">")) {
      # line starts a new record
      if (!is.null(header)) {
        # add the current record if it exists
        data$raw_seq[[length(data$raw_seq)+1]] <- list(sequence = sequence, header = header)
      }
      # start next record
      header <- sub("^>", "", line)
      sequence <- ""
    }
    else {
      sequence <- paste0(sequence, line)
    }
  }
  close(con)
  
  # handle last entry
  data$raw_seq[[length(data$raw_seq)+1]] <- list(sequence = sequence, header = header)
  
}

parse.fastq <- function(data) {
  # @arg data: Data object containing absolute or relative path to a FASTQ file
  # populates Data object's raw_seq list with (header, sequence, quality scores) lists
  
  con <- file(data$path, "r")

  # prepare containers
  header <- NULL
  sequence <- ""
  quality <- vector()
  ln <- 0

  while( TRUE ) {
    line = readLines(con, n = 1)
    if (length(line) == 0) {
      break  # reached end of file
    }
    position <- ln %% 4
    
    if (position == 0 && startsWith(line, "@")) {
      if (!is.null(header)) {
        data$raw_seq[[length(data$raw_seq)+1]] <- list(sequence = sequence, header = header, quality = quality)
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
    else if (position == 3) {
      quality <- convert_quality(line)
    }
    else {
      stop(paste("ERROR: Failed to parse FASTQ at line:\n", line))
    }
    ln <- ln + 1
  }
  close(con)
  
  # handle last entry
  data$raw_seq[[length(data$raw_seq)+1]] <- list(sequence = sequence, header = header, quality = quality)
}

convert_quality <- function(line) {
  # TODO: accommodate other quality score conversions (e.g., Solexa)
  # @arg line: string of encoded quality scores
  # @return vector of converted integer values
  
  # remove trailing whitespace (incl. line break)
  # split string into character vector
  line <- strsplit(sub("\\s+$", "", line), "")[[1]]
  
  result <- vector() 
  
  for (letter in line) {
    score <- utf8ToInt(letter) - 33
    if (score < 0 | score > 41) {
      stop(paste("ERROR: Unexpected integer value in convert_quality():", score))
    } else {
      result <- c(result, score)
    }
  }
  result
}

parse.default <- function(data) {
  warning(paste("highlineR does not know how to handle files of type ",
                class(data),
                "and can only be used on fasta and fastq files"))
}
