require(RUnit, quietly=TRUE)
source("import.R")

test.import <- function() {
  base <- "/home/lisamonique/Documents/highlineR/tests/test_data/"
  test <- paste0(base, "test.fa")
  
  # non existent file
  checkException(import(paste0(base, "test.fq")), "Error: file /home/lisamonique/Documents/highlineR/data/test.fq not found")
  # file with invalid extension
  checkException(import(paste0(base, "test.txt")), "Error: highlineR does not know how to handle files of type txt and can only be used on fasta and fastq files")
  # file with valid extension, invalid file type forced
  checkException(import(test, datatype = "txt"), "Error: highlineR does not know how to handle files of type txt and can only be used on fasta and fastq files")
  
  
  expected_Data <- structure(
    as.environment(
      list(
        path = test, 
        raw_seq = list(), 
        compressed = new.env()
      )
    ), 
    class = c("fasta", "Data", "environment")
  )
  checkEquals(Data(test), expected_Data)
  
  expected_session <- structure(new.env(), class = c("session", "environment"))
  
  # initialise empty session
  result <- init("test_sess")
  checkEquals(expected_session, result)
  
  assign(test, expected_Data, envir = expected_session)

  # import file to already initialised, empty session
  checkEquals(expected_session, import(test, session = "test_sess"))
  
  if (exists("test_sess")){
    rm(test_sess)
  }
  
  # import file to uninitalised session
  checkEquals(expected_session, import(test, session = "test_sess"))
  
  # import already imported file, expect warning, no change
  checkEquals(expected_session, import(test, session = "test_sess"))
  
  assign(paste0(base, "test2.fa"), Data(paste0(base, "test2.fa")), envir = expected_session)
  
  # import directory
  checkEquals(expected_session, import(base, session = "test_sess"))
}

test.parser <- function() {
  base <- "/home/lisamonique/Documents/highlineR/tests/test_data/"
  init("test_sess")
  expected_session <- import(base, session = "test_sess")
  
  # parse fasta file
  # parse(expected_session$`/home/lisamonique/Documents/highlineR/tests/test_data/test.fa`)
  # parse fastq file
  # parse already parsed file
  # parse entire session
}