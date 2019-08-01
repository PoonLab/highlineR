#' highlineR plotting
#' 
#' \code{highline} generates \href{https://www.hiv.lanl.gov/content/sequence/HIGHLIGHT/highlighter_top.html}{Highlighter}-style 
#' visualizations from one or more sequence alignments.  These plots use colour marks 
#' to annotate the differences of each sequence in the alignment from a reference sequence. 
#' 
#' @examples 
#'  highline(file, datatype='fasta', seqtype='nucleotide', unique=FALSE, 
#'  rf=1, mode='mismatch', sort_by='similarity', quiet=TRUE)
#' 
#' @param file a character vector containing one or more relative or absolute paths 
#' to the files containing sequence alignment(s).
#' @param datatype: specify file format - one of 'fasta', 'fastq' or 'csv'
#' @param seqtype: 'nucleotide' or 'amino acid'
#' @param rf: reading frame (defaults to 1) - affects 'svn' mode only
highline <- function(file, datatype='fasta', seqtype='nucleotide', unique=FALSE,
                     master=NA, rf=1, mode='mismatch', sort_by='similarity', 
                     use_sample=FALSE, quiet=TRUE) {

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
  plot(ses, mode=mode, rf=rf, sort_by=sort_by, use_sample=use_sample, 
       quiet=quiet)
}
