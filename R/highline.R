highline <- function(file, datatype='fasta', seqtype='nucleotide', unique=FALSE,
                     rf=1, mode='mismatch', sort_by='similarity', quiet=TRUE) {
  # A wrapper function to simplify the standard highlineR workflow
  #
  # Args:
  #   file: character vector containing one or more paths to alignment 
  #         files
  #   datatype: specify file format - one of 'fasta', 'fastq' or 'csv'
  #   seqtype: 'nucleotide' or 'amino acid'
  #   rf: reading frame (default 1) - affects 'tvt' mode only
  #   unique: if FALSE, gather identical sequences into unique variants
  #   mode: 'mismatch', 'svn' or 'tvt'
  #   sort_by: arrange sequences on plot by 'similarity' or 'frequency'
  #   quiet: if FALSE, suppress verbose output messages
  
  # input checks
  if (!is.character(file)) {
    stop(paste("Argument 'file' must be a character vector"))
  }
  if (any(!file.exists(file))) {
    stop(paste("File does not exist, stopping"))
  }
  
  # initialize highlineR session
  ses <- init_session()
  
  # import data from file
  for (f in file) {
    ses <- import_raw_seq(f, datatype=datatype, seqtype=seqtype, session=ses)
  }
  
  parse_raw_seq(ses)
  compress(ses, unique=unique)
  plot(ses, mode=mode, rf=rf, sort_by=sort_by, quiet=quiet)
}
