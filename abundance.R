compress <- function(x, ...){
  UseMethod("compress", x)
}


compress.session <- function(session, reads = 500, M = 10, unique = F, force = F) {
  # @arg session: environment containing imported and parsed sequence Data objects
  # @arg reads: expected number of reads
  # @arg M: number of random samples to average
  # populates each Data object's compressed environment with calculated sequence abundances
  
  # validate session contains imported sequences
  if (length(ls(session))==0) {
    stop("Error: session empty. Run highlineR::import(...)")
  }
  
  # calculate sequence abundance for each Data object
  for (data in ls(session)) {
    compress(get(data, envir = session, inherits = FALSE), M = M, reads = reads, unique = unique, force = force)
  }
}

compress.Data <- function(data, reads = 500, M = 10, unique = F, force = F) {
  # @arg data: Data object containing parsed list of sequences
  # @arg reads: expected number of reads
  # @arg M: number of random samples to average
  # populates Data object's compressed environment with (sequence: abunance) key:value pairs
  
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
    read_sample(data, M = M, reads = reads)
  }
}


compress.default <- function(x){
  warning(paste("Error: highlineR cannot compress objects of class ",
                class(x),
                "and can only be used on sessions or Data objects."))
}

read_sample <- function(data, reads, M) {
  # @arg data: Data object from which to sample
  # @arg reads: expected number of reads
  # @arg M: number of random samples to average
  # samples`reads` number of sequences from compressed sequences M times and stores average result in data$sample

  rm(list = ls(data$sample), envir = data$sample)
  dt <- data.frame(matrix(ncol = 3, nrow = length(data$compressed)))
  colnames(dt) <- c("seq", "reads", "f")
  dt$seq <- ls(data$compressed)
  for (i in 1:nrow(dt)) {
    dt[i, "reads"] <- data$compressed[[ls(data$compressed)[i]]]
  }
  N <- sum(dt$reads)
  dt$f <- dt$reads/N
  res <- data.frame(matrix(nrow = nrow(dt), ncol = M))
  rownames(res) <- dt$seq
  for (i in 1:M){
    samp <- sample(dt$seq, reads, prob = dt$f, replace = T)
    comp_samp <- as.data.frame(table(samp), stringsAsFactors = F)
    res[comp_samp$samp,i] <- comp_samp$Freq
  }
  ave <- round(rowMeans(res[complete.cases(res),], na.rm = T))
  list2env(as.list(ave), envir = data$sample)
}

resample <- function(x, ...) {
  UseMethod("resample", x)
}

resample.Data <- function(data, reads = 500, M = 10) {
  # @arg data: Data object from which to sample
  # @arg reads: expected number of reads
  # @arg M: number of random samples to average
  # resamples compressed sequences without re-compresing
  
  if (missing(reads) && missing(M)) {
    read_sample(data)
  }
  else if (missing(reads)) {
    read_sample(data, M = M)
  }
  else if (missing(M)) {
    read_sample(data, reads = reads)
  }
  else {
    read_sample(data, M = M, reads = reads)
  }
}

resample.session <- function(session, reads = 500, M = 10) {
  # @arg session: environment containing imported and parsed sequence Data objects
  # resamples all Data objects in session
  eapply(session, function(x) resample(x, reads = reads, M = M))
}
