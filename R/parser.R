#' Parse generic
#'
#' Parses the NGS file specified by the path variable of \code{Data} objects.
#'
#' @param x Data object or session object containing Data objects to parse.
#' @param encoding A character string representing the quality score encoding of FASTQ files. Options: "Sanger", "Solexa", "Illumina1.3", "Illumina1.5". "Illumina1.8".
#'
#' @return Parses Data objects by populating the object's \code{raw_seq} list with (header, sequence, [quality scores]) lists.
#' @export
parse_raw_seq <- function(x, ...) {
  UseMethod("parse_raw_seq", x)
}

#' @describeIn parse_raw_seq Parse all Data objects in session objects.
#'
#' @method parse_raw_seq session
#' @export
parse_raw_seq.session <- function(x, encoding = "sanger", ...) {
  for (data in ls(x)) {
    parse_raw_seq(get(data, envir = x, inherits = FALSE), encoding = encoding)
  }
}

#' @describeIn parse_raw_seq Parse Data objects of FASTA NGS files.
#'
#' @method parse_raw_seq fasta
#' @export
parse_raw_seq.fasta <- function(x, encoding = NULL, ...) {
  data <- x

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

    l <- -1 # sequence length

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
          if (nchar(sequence) > l){
            l <- nchar(sequence)
          }
        }
        # start next record
        header <- sub("^>", "", line)
        sequence <- ""
      }
      else if (inherits(data, "nucleotide") && grepl("^[ATGC-]", line, ignore.case = T)) {
        sequence <- paste0(sequence, toupper(line))
      }
      else if (inherits(data, "amino acid") && grepl("^[HMCDEILVAGSTKNQRFWYP-]", line, ignore.case = T)) {
        sequence <- paste0(sequence, toupper(line))
      }
      else if (line == "") {
        next
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

    # ensure sequences same length by adding gaps
    for (i in 1:(length(data$raw_seq))) {
      s <- data$raw_seq[[i]]
      if (nchar(s$sequence) < l) {
        data$raw_seq[[i]]$sequence <- paste0(s$sequence, paste(rep("-", l-nchar(s$sequence)), collapse = ""))
      }
    }
  }
}

#' @describeIn parse_raw_seq Parse Data objects of FASTQ NGS files.
#'
#' @method parse_raw_seq fastq
#' @export
parse_raw_seq.fastq <- function(x, encoding = "sanger", ...) {
  data <- x

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
            sequence <- toupper(line)
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

#' Convert encoded quality scores
#'
#' \code{convert_quality} converts encoded quality scores.
#'
#' @param line Character string of encoded quality scores.
#' @param encoding Character string of quality score encoding. Options: "Sanger" (default), "Solexa", "Illumina1.3", "Illumina1.5". "Illumina1.8".
#'
#' @return Numeric vector of converted quality scores.
#' @seealso \link{parse_raw_seq}
#' @export
convert_quality <- function(line, encoding = "sanger", ...) {
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

#' @describeIn parse_raw_seq Parse Data objects of custom, csv, NGS files.
#'
#' @method parse_raw_seq csv
#'
#' @keywords internal
parse_raw_seq.csv <- function(x, encoding = NULL, ...) {
  data <- x

  dt <- read.csv(data$path, sep=",", header = T, stringsAsFactors = F)
  l <- -1 # sequence length

  for (i in 1:nrow(dt)) {
    header <- paste(dt[i, "refname"], dt[i, "count"], sep = "_")
    # add gaps based on alignment offset
    sequence <- paste0(paste(rep("-", dt[i, "offset"] - (min(dt$offset))), collapse = ""), toupper(dt[i, "seq"]))

    # ignore sequences with ambiguous bases
    if (grepl("N", toupper(sequence))) {
      next
    }
    # keep track of maximum sequence length
    if (nchar(sequence) > l){
      l <- nchar(sequence)
    }
    data$raw_seq[[length(data$raw_seq)+1]] <- list(header = header, sequence = sequence)
  }
  # ensure sequences are the same length
  for (i in 1:(length(data$raw_seq))) {
    s <- data$raw_seq[[i]]
    if (nchar(s$sequence) < l) {
      data$raw_seq[[i]]$sequence <- paste0(s$sequence, paste(rep("-", l-nchar(s$sequence)), collapse = ""))
    }
  }
}

#' @describeIn parse_raw_seq default.
#'
#' @method parse_raw_seq default
#' @export
parse_raw_seq.default <- function(x) {
  warning(paste("highlineR does not know how to handle files of type ",
                class(x),
                "and can only be used on fasta and fastq files"))
}
