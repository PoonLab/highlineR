#' Compress generic
#'
#' Compresses the NGS file specified by the path variable of \code{Data} objects.
#'
#' @param x Data object or session object containing Data objects to compress.
#' @param N Integer expected number of reads.
#' @param M Integer number of random samples to select from dataset for read averaging.
#' @param unique A logical value specifying whether the NGS data file contains unique 
#'               sequences. If \code{True}, then it is assumed that read counts are 
#'               present in sequence labels. If \code{False}, then read counts are 
#'               inferred from the number of distinct instances of sequences in the file.
#' @param force A logical value specifying if previously imported file should be overwritten and re-imported.
#'
#' @return Compresses Data objects by populating the object's \code{compressed} environment with (sequence: abundance) key:value pairs.
#' @name compress
NULL

#' @rdname compress
#' @export
compress <- function(x, N, M, unique, force){
  UseMethod("compress", x)
}

#' @rdname compress
#' @export
compress.session <- function(x, N = 500, M = 10, unique = F, force = F) {
  # validate session contains imported sequences
  if (length(ls(x))==0) {
    stop("Error: session empty. Run highlineR::import(...)")
  }

  # calculate sequence abundance for each Data object
  for (data in ls(x)) {
    compress(get(data, envir = x, inherits = FALSE), M = M, N = N, unique = unique, force = force)
  }
}

#' @rdname compress
#' @export
compress.Data <- function(x, N = 500, M = 10, unique = F, force = F) {
  data <- x
  stopifnot(is.logical(unique))

  # ignore unparsed files
  if (length(data$raw_seq) == 0) {
    warning(paste("File", data$path, "ignored. Run highlineR::parse(...)"))
  }

  # ignore already compressed files
  else if(length(data$compressed) != 0 && force == F) {
    warning(paste("File", data$path, "ignored. Already compressed."))
  }
  else {
    master <- list("", -1)

    if (unique == F) {
      for (s in data$raw_seq){
        sequence <- s$sequence
        if (exists(sequence, envir = data$compressed)){
          # if sequence already in structure, increment count
          data$compressed[[sequence]] <-data$compressed[[sequence]] + 1
        }
        else{
          # otherwise, add sequence and initiate count
          data$compressed[[sequence]] <- 1
        }
        if (data$compressed[[sequence]] > master[[2]]) {
          # identify most abundant sequence
          master <- list(sequence, data$compressed[[sequence]])
        }
      }
    }
    else {
      for (s in data$raw_seq) {
        header <- strsplit(s$header, "[_-]")[[1]]
        count <- strtoi(trimws(header[length(header)]))
        sequence <- s$sequence

        if (exists(sequence, envir = data$compressed)){
          # if sequence already in structure, increment count
          data$compressed[[sequence]] <- data$compressed[[sequence]] + count
        }
        else{
          # otherwise, add sequence and initiate count
          data$compressed[[sequence]] <- count
        }

        if (data$compressed[[sequence]] > master[[2]]) {
          # identify most abundant sequence
          master <- list(sequence, data$compressed[[sequence]])
        }

      }
    }


    data$master <- master[[1]]
    read_sample(data, M = M, N = N)
  }
}

#' @rdname compress
#' @export
compress.default <- function(x, N, M, unique, force){
  warning(paste("Error: highlineR cannot compress objects of class ",
                class(x),
                "and can only be used on sessions or Data objects."))
}

#' @description \code{read_sample} samples \code{N} number of sequences from the compressed environment M times.
#' @return \code{read_sample} samples \code{N} number of sequences from the compressed evironment M times and stores the average result in the Data object's \code{sample} environment.
#' @rdname compress
#' @export
read_sample <- function(x, N, M) {
  data <- x

  rm(list = ls(data$sample), envir = data$sample)

  # calculate observed frequencies of variants
  dt <- data.frame(matrix(ncol = 3, nrow = length(data$compressed)))
  colnames(dt) <- c("seq", "N", "f")
  dt$seq <- ls(data$compressed)
  for (i in 1:nrow(dt)) {
    dt[i, "N"] <- data$compressed[[ls(data$compressed)[i]]]
  }
  total <- sum(dt$N)
  dt$f <- dt$N/total

  # create dataframe of sequenes * observed frequencies in M replicates
  res <- data.frame(matrix(nrow = nrow(dt), ncol = M))
  rownames(res) <- dt$seq
  for (i in 1:M){
    samp <- sample(dt$seq, N, prob = dt$f, replace = T)
    comp_samp <- as.data.frame(table(samp), stringsAsFactors = F)
    res[comp_samp$samp,i] <- comp_samp$Freq
  }
  # calculate average frequency of variants randomly selected in at least 75% of replicates
  ave <- round(rowMeans(res[rowSums(is.na(res)) < (0.25 * ncol(res)),], na.rm = T))
  s <- list2env(as.list(ave), envir = data$sample)

  s
}

#' @description \code{resample} resamples Data objects without re-compressing.
#' @return \code{resample} resamples Data objects by overwriting the object's \code{sample} environment with (sequence: abundance) key:value pairs without re-compressing.
#' @rdname compress
#' @export
resample <- function(x, N, M) {
  UseMethod("resample", x)
}

#' @rdname compress
#' @export
resample.session <- function(x, N = 500, M = 10) {
  eapply(x, function(x) resample(x, N = N, M = M))
}

#' @rdname compress
#' @export
resample.Data <- function(x, N = 500, M = 10) {
  data <- x

  if (missing(N) && missing(M)) {
    read_sample(data)
  }
  else if (missing(N)) {
    read_sample(data, M = M)
  }
  else if (missing(M)) {
    read_sample(data, N = N)
  }
  else {
    read_sample(data, M = M, N = N)
  }
}

