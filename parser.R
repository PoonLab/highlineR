#TODO: verify parsed list empty before parsing
#TODO: handle other type of quality scoring

parse <- function(x, ...) {
  UseMethod("parse", x)
}

parse.session <- function(session, quality_scoring = NULL, ...) {
  #arg session: environment containing imported sequence Data objects
  
  for (data in ls(session)) {
    parse(get(data, envir = session, inherits = FALSE), quality_scoring = quality_scoring)
  }
}

parse.fasta <- function(data, ...) {
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

parse.fastq <- function(data, quality_scoring = "sanger", ...) {
  # @arg data: Data object containing absolute or relative path to a FASTQ file
  # populates Data object's raw_seq list with (header, sequence, quality scores) lists
  
  #validate quality scoring method
  quality_scoring <- match.arg(tolower(quality_scoring), c("sanger", "solexa"))
  
  con <- file(data$path, "r")

  # prepare containers
  header <- NULL
  sequence <- ""
  quality <- vector()
  ln <- 0
  res <- NULL
  
  res <- tryCatch({
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
        quality <- convert_quality(line, quality_scoring = quality_scoring)
      }
      else {
        stop(paste("ERROR: Failed to parse FASTQ at line:\n", line))
      }
      ln <- ln + 1
    }
    # handle last entry
    data$raw_seq[[length(data$raw_seq)+1]] <- list(sequence = sequence, header = header, quality = quality)
  },
  error = function(e) {
    warning(e)
    data$raw_seq <- list()
  },
  warning = function(w) {
    warning(w)
  },
  finally = {
    close(con)
  })
  if(inherits(res, "error")){
    print(res)
  }
}

convert_quality <- function(line, quality_scoring = "sanger", ...) {
  # @arg line: string of encoded quality scores
  # @return vector of converted integer values
  
  # remove trailing whitespace (incl. line break)
  # split string into character vector
  line <- strsplit(sub("\\s+$", "", line), "")[[1]]
  
  result <- vector() 
  
  
  for (letter in line) {
    if (quality_scoring == "sanger") {
      score <- utf8ToInt(letter) - 33
      min <- 0
    }
    else if(quality_scoring == "solexa") {
      score <- utf8ToInt(letter) - 59
      min <- -5
    }
    
    if (score < min | score > 41) {
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
