highline <- function(fn, ..., datatype='fasta', seqtype='nucleotide') {
  # A wrapper function to simplify the standard workflow
  
  # initialize highlineR session
  ses <- init_session()
  
  # import data from file
  ses <- import_raw_seq(fn, datatype=datatype, seqtype=seqtype, session=ses)
  
  parse_raw_seq(ses)
  compress(ses, unique=FALSE)
  plot(ses)
}

fn <- '~/git/highlineR/working/04013448.fasta'
