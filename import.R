# TODO: nucleotide or amino acid data type
# TODO: import_raw_seq(vector)

Data <-  function(path, datatype = tail(strsplit(path, "\\.")[[1]], n = 1)) {
  # @arg path absolute path to sequence file
  # @arg datatype file type, default: file extension
  # @return s3 Data object to hold raw and processed sequencing data for single file
  
  # validate file exists
  if (!file.exists(path)) {
    stop(paste0("ERROR: file '", path, "' not found")
    )
  }
  
  # validate file type
  if (tolower(datatype) %in% c("fasta", "fa")) {
    datatype <- "fasta"
  }
  else if (tolower(datatype) %in% c("fastq", "fq")) {
    datatype <- "fastq"
  }
  else {
    stop(paste0("ERROR: file '", path, "' not imported. highlineR does not know how to handle files of type ",
               datatype,
               " and can only be used on fasta and fastq files")
    )
  }
  
  # create s3 Data structure
  data <- structure(
    as.environment(
      list(
        path = path, 
        raw_seq = list(), 
        compressed = 
          structure(
            new.env(),
            class = c("compressed", "environment")
            ),
        master = character(),
        seq_diff = matrix()
        )
      ), 
    class = c(datatype, "Data", "environment")
    )
  
  parent.env(data$compressed) <- data
  
  data
}

remove_Data <- function(data, session = "highlineR.session") {
  # @arg data Data object or string name of Data object to be removed
  # @arg session session containing Data object to be removed
  # removes specified Data object from specified session environment
  if (is.character(session)){
    session <- get(session)
  }
  
  if (is.character(data)) {
    # check if specified Data objec exists in specified session
    stopifnot(exists(data, envir = session))
    # remove
    rm(list=data, envir = session)
  }
  else {
    # if Data object provided
    stopifnot(inherits(data, "Data"))
    rm(list=data$path, envir = parent.env(data))
  }
}


init_session <- function(session = "highlineR.session") {
  # @arg session string name for environment
  # @return empty environment to hold sequence data
  
  assign(session, 
         structure(new.env(), class = c("session", "environment")),
         envir = .GlobalEnv
  )
  get(session)
}

close_session <- function(session = "highlineR.session") {
  # @arg session string name or session object
  # removes session environment object
  
  if (! is.character(session)) {
    session <- deparse(substitute(session))
  }
  stopifnot(exists(session))
  
  # check session argument is of session class
  stopifnot(inherits(get(session, envir = .GlobalEnv), "session"))
  
  rm(list=paste(session), envir = .GlobalEnv)
}

import_file <- function(path, datatype, session, force = FALSE) {
  # @arg path: absolute path to sequence file
  # @arg datatype: file type, optional
  # @arg session: string name of environment to load sequence files into
  # imports file into specified session
  
  if (missing(session)) {
    stop("highlineR session not specified.")
  }
  stopifnot(is.character(session))
  
  stopifnot(is.character(path))
  stopifnot(file.exists(path))
  
  if (! missing(datatype)) {
    stopifnot(is.character(datatype))
  }
  
  # create environment if it doesn't exist
  if (! exists(session)) {
    init_session(session)
  }
  
  # ignore files already imported unless forced
  if (force == FALSE && exists(path, envir = get(session), inherits = FALSE)) {
    warning(paste("File", path, "ignored. Already imported in", session, "session."))
  }
  # otherwise create Data object within specified session
  else {
    if (missing(datatype)) {
      data <- Data(path)
    }
    else {
      data <- Data(path, datatype = datatype)
    }
    parent.env(data) <- get(session)
    assign(path, data, envir = get(session))
  }
}


import_raw_seq <- function(path, datatype, session = "highlineR.session", force = FALSE) {
  # @arg path: absolute path to sequence file or directory containing sequence files
  # @arg datatype: file type, optional
  # @arg session: string name of environment to load sequence files to, default highlineR.data
  # @return environment containing imported sequence file(s) 

  # validate path
  if(!file.exists(path)) {
    stop(paste0("ERROR: path '", path, "' not valid"))
  }
  # stopifnot(file.exists(path))
  
  # if path is directory, import each file
  if (dir.exists(path)) {
    # remove trailing backslash
    path <- sub("\\/$", "", path)
    tmp.list.1 <- list.files(path, full.names = TRUE, recursive = TRUE)
    
    for (i in 1:length(tmp.list.1)) {
      result <- tryCatch(
        if (missing(datatype)) {
          import_file(tmp.list.1[i], session = session, force = force)
        }
        else {
          import_file(tmp.list.1[i], datatype = datatype, session = session, force = force)
        }
      , error = function(e) {
        warning(e)
        e
      }
      , warning = function(w) {
        warning(w)
      }
      )
      if (inherits(result, "error"))
        next
    }
  }
  # if path is file, import file
  else {
    if (missing(datatype)) {
      import_file(path, session = session, force = force)
    }
    else {
      import_file(path, datatype = datatype, session = session, force = force)
    }
    
  }
  
  get(session)
}
