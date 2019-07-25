highline <- function(file, datatype='fasta', seqtype='nucleotide',
                     mode='mismatch', sort_by='similarity', quiet=TRUE) {
  # A wrapper function to simplify the standard workflow
  #
  # Args:
  #   file: character vector containing one or more paths to alignment 
  #         files
  #   datatype: specify file format - one of 'fasta', 'fastq' or 'csv'
  #   seqtype: 'nucleotide' or 'amino acid'
  #   mode: 'mismatch', 'svn' or 'tvt'
  #   sort_by: arrange sequences on plot by 'similarity' or 'frequency'
  #   quiet: if FALSE, suppress verbose output messages
  
  # input checks
  if (!is.character(file)) {
    stop(paste("Argument 'file' must be a character vector"))
  }
  if (any(!file.exists(file))) {
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
