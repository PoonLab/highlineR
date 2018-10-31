new_Data <-  function(file, datatype=tail(strsplit(file, "\\.")[[1]], n=1)){
  #@arg file absolute path to sequence file
  #@arg datatype file type
  #@return data s3 object to hold raw and processed sequencing data for single plot
  
  #validate file exists
  if (!file.exists(file)){
    stop(paste("Error: file", file, "not found"),
         call. = FALSE
    )
  }
  
  #validate file type
  if (tolower(datatype) %in% c("fasta", "fa")){
    datatype = "fasta"
  }
  else if (tolower(datatype) %in% c("fastq", "fq")){
    datatype = "fastq"
  }
  else {
    stop(paste("highlineR does not know how to handle files of type ",
               datatype,
               "and can only be used on fasta and fastq files"),
         call. = FALSE
    )
  }
  
  structure(list(path = file, raw_seq = list(), compressed = new.env()), 
            class = c(datatype, "Data")
  )
}

init <- function(session = character()){
  #initiate empty environment to hold sequence data
  #@arg session string name for environment
  
  assign(session, new.env(), envir = .GlobalEnv)
}

import <- function(path, session="highlineR.data"){
  #@arg path absolute path to sequence file or directory containing sequence files
  #@arg session name of environment to load sequence files to, default highlineR.data
  #@return none
  
  #create environment if doesn't exist
  if (! exists(session)){
    init(session)
  }
  
  #validate path
  stopifnot(file.exists(path))
  
  #if path is directory, import each file
  if (dir.exists(path)){
    tmp.list.1 <- list.files(path, full.names = TRUE)
    for (i in 1:length(tmp.list.1)){
      import_file(tmp.list.1[i], session=session)
    }
  }
  else {
    import_file(path, session=session)
  }
}

import_file <- function(path, session){
  #@arg path absolute path to sequence file or directory containing sequence files
  #@arg session name of environment to load sequence files to, default highlineR.data
  #@return none
  
  #ignore files already imported
  if (exists(path, envir = get(session))){
    warning(paste("File", path, "ignored. Already imported in", session, "session."))
  }
  else{
    # otherwise create data object
    data <- new_Data(path)
    assign(path, data, envir = get(session))
  }
}
