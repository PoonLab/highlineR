highline <- function(file, ..., datatype='fasta', seqtype='nucleotide',
                     mode='mismatch', sort_by=quiet=FALSE) {
  # A wrapper function to simplify the standard workflow
  #
  # Args:
  #   file: path to alignment file (format FASTA or FASTQ, see <datatype>)
  #         or a character vector of multiple paths
  #   ...: additional file paths
  #   
  
  # input checks
  if (!file_test(file)) {
    stop(paste("File ", file, " does not exist, stopping"))
  }
  
  # initialize highlineR session
  ses <- init_session()
  
  # import data from file
  if (length(fn) == 1 && )
  files <- c(list(fn), list(...))
  for (file in files) {
    ses <- import_raw_seq(file, datatype=datatype, seqtype=seqtype, session=ses)
  }
  
  parse_raw_seq(ses)
  compress(ses, unique=FALSE)
  plot(ses, quiet=quiet)
}

fn <- '~/git/highlineR/working/04013448.fasta'

files <- Sys.glob('~/git/highlineR/working/*.fasta')
