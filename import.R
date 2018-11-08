#TODO: is.Data, is.fasta, is.fastq etc
#TODO: close(sesion)

Data <-  function(file, datatype = tail(strsplit(file, "\\.")[[1]], n = 1)) {
  # @arg file absolute path to sequence file
  # @arg datatype file type, default: file extension
  # @return s3 Data object to hold raw and processed sequencing data for single file
  
  # validate file exists
  if (!file.exists(file)) {
    stop(paste("Error: file", file, "not found"),
         call. = FALSE
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
    stop(paste("Error: file", file, "not imported. highlineR does not know how to handle files of type",
               datatype,
               "and can only be used on fasta and fastq files"),
         call. = FALSE
    )
  }
  
  # create s3 Data structure
  data <- structure(
    as.environment(
      list(
        path = file, 
        raw_seq = list(), 
        compressed = new.env(),
        master = character(),
        seq_diff = new.env()
        )
      ), 
    class = c(datatype, "Data", "environment")
    )
  
  
  data
}


init <- function(session = "highlineR.data") {
  # @arg session string name for environment
  # @return empty environment to hold sequence data
  
  assign(session, 
         structure(new.env(), class = c("session", "environment")),
         envir = .GlobalEnv
  )
  get(session)
}


import_file <- function(file, datatype, session) {
  # @arg file: absolute path to sequence file
  # @arg datatype: file type, optional
  # @arg session: string name of environment to load sequence files into
  # imports file into specified session
  
  stopifnot(is.character(file))
  if (! missing(datatype)) {
    stopifnot(is.character(datatype))
  }
  if(! missing(session)) {
    stopifnot(is.character(session))
  }
  stopifnot(file.exists(file))
  
  # create environment if it doesn't exist
  if (! exists(session)) {
    init(session)
  }
  
  # ignore files already imported
  if (exists(file, envir = get(session), inherits = FALSE)) {
    warning(paste("File", file, "ignored. Already imported in", session, "session."))
  }
  # otherwise create Data object within specified session
  else {
    if (missing(datatype)) {
      assign(file, Data(file), envir = get(session))
    }
    else {
      assign(file, Data(file, datatype = datatype), envir = get(session))
    }
  }
}


import <- function(path, datatype, session = "highlineR.data") {
  # @arg path: absolute path to sequence file or directory containing sequence files
  # @arg datatype: file type, optional
  # @arg session: string name of environment to load sequence files to, default highlineR.data
  # @return environment containing imported sequence file(s) 

  # validate path
  stopifnot(file.exists(path))
  
  # if path is directory, import each file
  if (dir.exists(path)) {
    # remove trailing backslash
    path <- sub("\\/$", "", path)
    tmp.list.1 <- list.files(path, full.names = TRUE)
    
    for (i in 1:length(tmp.list.1)) {
      result <- tryCatch(
        if (missing(datatype)) {
          import_file(tmp.list.1[i], session = session)
        }
        else {
          import_file(tmp.list.1[i], datatype = datatype, session = session)
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
      import_file(path, session = session)
    }
    else {
      import_file(path, datatype = datatype, session = session)
    }
    
  }
  
  get(session)
}
