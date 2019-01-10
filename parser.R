parse <- function(x, ...) {
  UseMethod("parse", x)
}

parse.session <- function(session, encoding = "sanger", ...) {
  #arg session: environment containing imported sequence Data objects
  
  for (data in ls(session)) {
    parse(get(data, envir = session, inherits = FALSE), encoding = encoding)
  }
}

parse.fasta <- function(data, encoding = NULL, ...) {
  # @arg data: Data object containing absolute or relative path to a FASTA file
  # populates Data object's raw_seq list with (header, sequence) lists
  
  # ignore already parsed files
  if (length(data$raw_seq) != 0) {
    warning(paste0("ERROR: file '", data$path, "' ignored. Already parsed."))
  }
  else{
    con <- file(data$path, "r")
    
    # prepare containers
    header <- NULL
    sequence <- ""
    ln <- 0
    
    while( TRUE ) {
      line = readLines(con, n = 1)
      if (length(line) == 0) {
        break  # reached end of file
      }
      
      if (startsWith(line, ">")) {
        # line starts a new record
        if (!is.null(header)) {
          # add the current record if it exists
          data$raw_seq[[length(data$raw_seq)+1]] <- list(header = header, sequence = sequence)
        }
        # start next record
        header <- sub("^>", "", line)
        sequence <- ""
      }
      else if (inherits(data, "nucleotide") && grepl("^[A|T|G|C|-]", line)) {
        sequence <- paste0(sequence, line)
      }
      else if (inherits(data, "amino acid") && grepl("^[H|M|C|D|E|I|L|V|A|G|S|T|K|N|Q|R|F|W|Y|P|-]", line)) {
        sequence <- paste0(sequence, line)
      }
      else {
        # incorrect format
        stop(paste("ERROR: Failed to parse FASTA at line:", ln))
      }
      ln <- ln + 1
    }
    close(con)
    
    # handle last entry
    data$raw_seq[[length(data$raw_seq)+1]] <- list(header = header, sequence = sequence)
  }
}

parse.fastq <- function(data, encoding = "sanger", ...) {
  # @arg data: Data object containing absolute or relative path to a FASTQ file
  # @arg encoding FASTQ file quality score encoding. Options: "Sanger", "Solexa", "Illumina1.3", "Illumina1.5". "Illumina1.8"
  # populates Data object's raw_seq list with (header, sequence, quality scores) lists
  
  # ignore already parsed files
  if (length(data$raw_seq) != 0) {
    warning(paste0("ERROR: file '", data$path, "' ignored. Already parsed."))
  }
  else {
    con <- file(data$path, "r")
    
    # prepare containers
    header <- NULL
    sequence <- ""
    quality <- vector()
    ln <- 0
    
    res <- tryCatch(
      {
        while( TRUE ) {
          line = readLines(con, n = 1)
          if (length(line) == 0) {
            break  # reached end of file
          }
          position <- ln %% 4
          
          if (position == 0 && startsWith(line, "@")) {
            if (!is.null(header)) {
              # add the current record if it exists
              data$raw_seq[[length(data$raw_seq)+1]] <- list(header = header, sequence = sequence, quality = quality)
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
            quality <- convert_quality(line, encoding = encoding)
          }
          else {
            # incorrect format
            stop(paste("ERROR: Failed to parse FASTQ at line:", ln))
          }
          ln <- ln + 1
        }
        # handle last entry
        data$raw_seq[[length(data$raw_seq)+1]] <- list(header = header, sequence = sequence, quality = quality)
      },
      error = function(e) {
        # if error in parsing, reset raw_seq list
        warning(e)
        data$raw_seq <- list()
      },
      warning = function(w) {
        warning(w)
      }
    )
    close(con)
    
    if(inherits(res, "error")){
      print(res)
    }
  }
}

convert_quality <- function(line, encoding = "sanger", ...) {
  # @arg line: string of encoded quality scores
  # @return vector of converted integer values
  
  #validate quality scoring method
  encoding <- match.arg(tolower(encoding), c("sanger", "solexa", "illumina1.3", "illumina1.5", "illumina1.8"))
  
  # remove trailing whitespace (incl. line break)
  # split string into character vector
  line <- strsplit(sub("\\s+$", "", line), "")[[1]]
  
  result <- vector() 
  
  min <- 0
  max <- 40
  for (letter in line) {
    if (encoding == "sanger") {
      score <- utf8ToInt(letter) - 33
    }
    else if(encoding == "solexa") {
      score <- utf8ToInt(letter) - 64
      min <- -5
    }
    else if(encoding == "illumina1.3") {
      score <- utf8ToInt(letter) - 64
    }
    else if(encoding == "illumina1.5") {
      score <- utf8ToInt(letter) - 64
      min <- 3
    }
    else if(encoding == "illumina1.8") {
      score <- utf8ToInt(letter) - 33
      max <- 41
    }
    
    if (score < min | score > max) {
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
