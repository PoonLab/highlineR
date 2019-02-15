compress <- function(x){
  UseMethod("compress", x)
}


compress.session <- function(session) {
  # @arg session: environment containing imported and parsed sequence Data objects
  # populates each Data object's compressed environment with calculated sequence abundances
  
  # validate session contains imported sequences
  if (length(ls(session))==0) {
    stop("Error: session empty. Run highlineR::import(...)")
  }
  
  # calculate sequence abundance for each Data object
  for (data in ls(session)) {
    compress(get(data, envir = session, inherits = FALSE))
  }
}


compress.Data <- function(data) {
  # @arg data: Data object containing parsed list of sequences
  # populates Data object's compressed environment with (sequence: abunance) key:value pairs
  
  # ignore unparsed files
  if (length(data$raw_seq) == 0) {
    warning(paste("File", data$path, "ignored. Run highlineR::parse(...)"))
  }
  # ignore already compressed files
  else if(length(data$compressed) != 0) {
    warning(paste("File", data$path, "ignored. Already compressed."))
  }
  else {
    master <- list("", -1)
    
    # for (s in data$raw_seq){
    #   sequence <- s$sequence
    #   if (exists(sequence, envir = data$compressed)){
    #     # if sequence already in structure, increment count
    #     data$compressed[[sequence]] <-data$compressed[[sequence]] + 1
    #   }
    #   else{
    #     # otherwise, add sequence and initiate count
    #     data$compressed[[sequence]] <- 0
    #   }
    #   if (data$compressed[[sequence]] > master[[2]]) {
    #     # identify most abundant sequence
    #     master <- list(sequence, data$compressed[[sequence]])
    #   }
    # }
    
    
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
    
    data$master <- master[[1]]
    read_sample(data)
  }
}


compress.default <- function(x){
  warning(paste("Error: highlineR cannot compress objects of class ",
                class(x),
                "and can only be used on sessions or Data objects."))
}

read_sample <- function(data) {
  rm(list = ls(data$sample), envir = data$sample)
  dt <- data.frame(matrix(ncol = 3, nrow = length(data$compressed)))
  colnames(dt) <- c("seq", "reads", "f")
  dt$seq <- ls(data$compressed)
  for (i in 1:nrow(dt)) {
    dt[i, "reads"] <- data$compressed[[ls(data$compressed)[i]]]
  }
  N <- sum(dt$reads)
  dt$f <- dt$reads/N
  samp <- sample(dt$seq, 10000, prob = dt$f, replace = T)
  list2env(as.list(table(samp)), envir = data$sample)
}

resample <- function(x) {
  UseMethod("resample", x)
}
resample.Data <- function(data) {
  read_sample(data)
}
resample.session <- function(session) {
  eapply(session, resample)
}
