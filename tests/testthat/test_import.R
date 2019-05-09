context("Importing")
base <- paste0(getwd(), "/test_data/")


test_that("invalid files do not open", {
  if (exists("test_sess", envir = .GlobalEnv)){
    rm(test_sess, envir = .GlobalEnv)
  }

  # non-existent file
  expect_error(import_raw_seq(paste0(base, "test.fq"), session = "test_sess"),
               paste0("ERROR: path '", paste0(base, "test.fq", "' not valid")))
  # file with invalid extension, default import
  expect_error(import_raw_seq(paste0(base, "invalid/test.txt"), session = "test_sess"),
               paste0("ERROR: file '", paste0(base, "invalid/test.txt", "' not imported. highlineR does not know how to handle files of type txt and can only be used on fasta and fastq files")))
  # file with valid extension, specify import with invalid extension
  expect_error(import_raw_seq(paste0(base, "valid/fasta/test.fa"), datatype = "txt", session = "test_sess"),
               paste0("ERROR: file '", paste0(base, "valid/fasta/test.fa", "' not imported. highlineR does not know how to handle files of type txt and can only be used on fasta and fastq files")))

})


test_that("highlineR s3 Data objects created correctly", {
  expected_Data <- structure(
    as.environment(
      list(
        path = paste0(base, "valid/fasta/test.fa"),
        raw_seq = list(),
        compressed =
          structure(
            new.env(),
            class = c("compressed", "environment")
          ),
        sample =
          structure(
            new.env(),
            class = c("compressed", "environment")
          ),
        master = character(),
        seq_diff = matrix()
      )
    ),
    class = c("fasta", "Data", "environment")
  )
  parent.env(expected_Data$compressed) <- expected_Data

  # create Data object
  expect_error(result <- Data(paste0(base, "valid/fasta/test.fa")),
               NA)
  # check correct structure with correct class
  expect_equal(expected_Data, result)
  expect_s3_class(result, "Data")
})


test_that("sessions initate correctly with or without file import", {
  test_filename <- paste0(base, "valid/fasta/test.fa")
  if (exists("test_sess", envir = .GlobalEnv)){
    rm(test_sess, envir = .GlobalEnv)
  }
  if (exists("test_sess_2", envir = .GlobalEnv)){
    rm(test_sess_2, envir = .GlobalEnv)
  }

  # empty session structure
  expected_session <- structure(new.env(), class = c("session", "environment"))

  # initialise empty session
  test_sess <- init_session()
  expect_equal(expected_session, test_sess)

  # session structure with single file
  expected_Data <- structure(
    as.environment(
      list(
        path = test_filename,
        raw_seq = list(),
        compressed =
          structure(
            new.env(),
            class = c("compressed", "environment")
          ),
        sample =
          structure(
            new.env(),
            class = c("compressed", "environment")
          ),
        master = character(),
        seq_diff = matrix()
      )
    ),
    class = c("fasta", "Data", "environment")
  )
  parent.env(expected_Data$compressed) <- expected_Data
  assign(test_filename, expected_Data, envir = expected_session)

  # import file to already initialised, empty session
  expect_equal(expected_session, import_raw_seq(test_filename, session = "test_sess"))
  if (exists("test_sess", envir = .GlobalEnv)){
    rm(test_sess, envir = .GlobalEnv)
  }

  # import file to uninitalised session
  test_sess_2 <<- import_raw_seq(test_filename)
  expect_equal(expected_session, test_sess_2)

  # import already imported file
  expect_warning(import_raw_seq(test_filename, session = "test_sess_2"),
                 paste0("File ", base, "valid/fasta/test.fa ignored. Already imported in session"))
  expect_equal(expected_session, test_sess_2) # no new file added

  # import already imported file, forced
  expect_warning(import_raw_seq(test_filename, session = "test_sess_2", force = TRUE),
                 NA)
  if (exists("test_sess_2", envir = .GlobalEnv)){
    rm(test_sess_2, envir = .GlobalEnv)
  }

  assign(paste0(base, "valid/fasta/test2.fa"), Data(paste0(base, "valid/fasta/test2.fa")), envir = expected_session)

  # import directory
  expect_error(test_sess_2 <- import_raw_seq(paste0(base, "valid/fasta/")),
                 NA)
  expect_equal(expected_session, test_sess_2)

  if (exists("test_sess_2", envir = .GlobalEnv)){
    rm(test_sess_2, envir = .GlobalEnv)
  }
})

test_that("sessions can be initialised and closed", {
  if (exists("test_sess", envir = .GlobalEnv)){
    rm(test_sess, envir = .GlobalEnv)
  }

  # initiate session
  expect_error(test_sess <<- init_session(),
               NA)
  # check if session exists and is of class session
  expect_true(exists("test_sess", envir = .GlobalEnv))
  expect_s3_class(test_sess, "session")

  # close session
  expect_error(close_session(test_sess),
               NA)
  # check that session no longer exists
  expect_false(exists("test_sess", envir = .GlobalEnv))

  if (exists("test_sess", envir = .GlobalEnv)){
    rm(test_sess, envir = .GlobalEnv)
  }
})

test_that("Data objects can be removed from sessions", {
  if (exists("test_sess", envir = .GlobalEnv)){
    rm(test_sess, envir = .GlobalEnv)
  }
  test_filename <- paste0(base, "valid/fasta/test.fa")

  # initiate session with test file
  test_sess <<- import_raw_seq(test_filename)

  # check that file exists in session
  expect_true(exists(test_filename, envir = test_sess))

  # remove file from session
  expect_error(remove_Data(test_filename, session = "test_sess"),
               NA)
  # check that file no longer exists in session
  expect_false(exists(test_filename, envir = test_sess))

  if (exists("test_sess", envir = .GlobalEnv)){
    rm(test_sess, envir = .GlobalEnv)
  }
})
