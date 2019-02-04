require(testthat, quietly=TRUE)
source("../import.R")
source("../parser.R")

test_that("fasta files parsed correctly", {
  base <- paste0(getwd(), "/test_data/valid/fasta/")
  test_filename <- paste0(base, "test.fa")
  
  # parse valid fasta data object
  test_Data <- Data(test_filename)
  parse(test_Data)
  
  expected_raw_seq <- lapply(
    list(
      list("adam", "ACTGACTAGCTAGCTAACTG"),
      list("bob", "GCATCGTAGCTAGCTACGAT"),
      list("courtney", "CATCGATCGTACGTACGTAG"),
      list("danny", "ATCGATCGATCGTACGATCG"),
      list("john", "ACTGACTAGCTAGCTAACTG"),
      list("sam", "ACTGACTAGCTAGCTAACTG")
    ),
    setNames,
    c("header", "sequence")
    )
  
  expect_equal(test_Data$raw_seq, expected_raw_seq)
  
  test_filename <- paste0(base, "test2.fa")
  # parse valid fasta data object with multi-line sequences
  test_Data <- Data(test_filename)
  expect_error(parse(test_Data),
               NA)
  expect_equal(test_Data$raw_seq, expected_raw_seq)
})

test_that("quality score conversion correct", {
  # default = sanger
  expect_equal(convert_quality(intToUtf8(33:73)),
               0:40)
  expect_equal(convert_quality(intToUtf8(33:73), encoding = "sanger"),
               0:40)
  expect_equal(convert_quality(intToUtf8(59:104), encoding = "solexa"),
               -5:40)
  expect_equal(convert_quality(intToUtf8(64:104), encoding = "illumina1.3"),
               0:40)
  expect_equal(convert_quality(intToUtf8(67:104), encoding = "illumina1.5"),
               3:40)
  expect_equal(convert_quality(intToUtf8(33:74), encoding = "illumina1.8"),
               0:41)
  # invalid quality scoring encoding
  expect_error(convert_quality(intToUtf8(33:73), encoding = "invalid"),
               "'arg' should be one of “sanger”, “solexa”, “illumina1.3”, “illumina1.5”, “illumina1.8”")
})

test_that("fastq files parsed correctly", {
  base <- paste0(getwd(), "/test_data/valid/fastq/")
  test_filename <- paste0(base, "test.fq")
  
  # parse valid fastq data object
  test_Data <- Data(test_filename)
  expect_error(parse(test_Data),
               NA)
  
  expected_raw_seq <- lapply(
    list(
      list("adam", "ACTGACTAGCTAGCTAACTG", convert_quality("9C;=;=<9@4868>9:67AA")),
      list("bob", "GCATCGTAGCTAGCTACGAT", convert_quality("1/04.72,(003,-2-22+0")),
      list("courtney", "CATCGATCGTACGTACGTAG", convert_quality("?7?AEEC@>=1?A?EEEB9E")),
      list("danny", "ATCGATCGATCGTACGATCG", convert_quality(">=2.660/?:36AD;0<147")),
      list("john", "CATCGATCGTACGTACGTAG", convert_quality("8;;;>DC@DAC=B?C@9?B?")),
      list("sam", "GCATCGTAGCTAGCTACGAT", convert_quality("-/CA:+<599803./2065?"))
    ),
    setNames,
    c("header", "sequence", "quality")
  )
  
  expect_equal(test_Data$raw_seq, expected_raw_seq)
})

test_that("invalid files not parsed", {
  base <- paste0(getwd(), "/test_data/valid/")
  
  # import fasta format file as fastq and parse
  test_filename <- paste0(base, "fasta/test.fa")
  test_Data <- Data(test_filename, datatype = "fq")
  expect_error(parse(test_Data), "ERROR: Failed to parse FASTQ at line: 0")
  
  # import fasta format file as fastq and parse
  test_filename <- paste0(base, "fastq/test.fq")
  test_Data <- Data(test_filename, datatype = "fa")
  expect_error(parse(test_Data), "ERROR: Failed to parse FASTA at line: 0")
  
  # try to parse already parsed file
  test_filename <- paste0(base, "fasta/test.fa")
  test_Data <- Data(test_filename)
  parse(test_Data)
  expect_warning(parse(test_Data), paste0("ERROR: file '", getwd(), "/test_data/valid/fasta/test.fa' ignored. Already parsed."))
})

test_that("directories can be parsed", {
  base <- paste0(getwd(), "/test_data/valid/fasta/")
  if (exists("test_sess")){
    rm (test_sess, envir = .GlobalEnv)
  }
  
  test <- import_raw_seq(base, session = "test_sess")
  expect_error(parse(test), NA)
  expect_setequal(ls(test), c(paste0(base, "test.fa"), paste0(base, "test2.fa")))
  
  if (exists("test_sess")){
    rm (test_sess, envir = .GlobalEnv)
  }
})
