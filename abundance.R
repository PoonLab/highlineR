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
    warning(paste("File", get(data, envir = session)$path, "ignored. Run highlineR::parse(...)"))
  }
  # ignore already compressed files
  else if(length(data$compressed) != 0) {
    warning(paste("File", get(data, envir = session)$path, "ignored. Already compressed."))
  }
  else {
    master <- list("", -1)
    
    for (s in data$raw_seq){
      sequence <- s$sequence
      if (exists(sequence, envir = data$compressed)){
        # if sequence already in structure, increment count
        data$compressed[[sequence]] <-data$compressed[[sequence]] + 1 
      }
      else{
        # otherwise, add sequence and initiate count
        data$compressed[[sequence]] <- 0
      }
      if (data$compressed[[sequence]] > master[[2]]) {
        # identify most abundant sequence
        master <- list(sequence, data$compressed[[sequence]])
      }
    }
    data$master <- master[[1]]
  }
}


compress.default <- function(x){
  warning(paste("Error: highlineR cannot compress objects of class ",
                class(x),
                "and can only be used on sessions or Data objects."))
}
