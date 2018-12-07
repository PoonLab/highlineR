summary.Data <- function(data, ...) {
  # @arg data highlineR Data object to be summarized
  # @return list of summarized variables
  
  fn <- data$path # filename
  dt <- class(data)[1] # datatype of file
  seq_num <- length(data$raw_seq) # number of sequences
  seq_len <- nchar(data$raw_seq[[1]]$sequence) # length of sequence
  var_num <- length(ls(data$compressed)) # number of variants
  var_counts <- sort(data$compressed) # variant counts sorted by abundance
  abun_var <- rownames(var_counts)[1] # most abundant variant
  
  res <- list(fn = fn,
              dt = dt,
              seq_num = seq_num,
              seq_len = seq_len,
              var_num = var_num,
              abun_var = abun_var,
              vars = var_counts)
  class(res) <- "summary.Data"
  res
}

print.summary.Data <- function(x, ...) {
  # print summary of Data object
  
  cat("File: ")
  cat(x$fn)
  
  cat("\nDatatype: ")
  cat(x$dt)
  
  cat("\nNumber of Sequences: ")
  cat(x$seq_num)
  
  cat("\nSequence Length: ")
  cat(x$seq_len)
  
  cat("\nNumber of Variants: ")
  cat(x$var_num)
  
  cat("\nMost Abundant Variant: ")
  cat(x$abun_var)
  cat("\n\nVariants:\n")
  
  print(x$vars)
  
  cat("\n")
  cat(rep("_", 40))
  cat("\n\n")
}

summary.session <- function(session, ...) {
  # @arg session highlineR session of sequence Data objects to be summarized
  # @return list of summarized variables for each sequence Data object
  
  res <- eapply(session, summary)
  class(res) <- "summary.session"
  res
}

print.summary.session <- function(x, ...) {
  # print summary of each Data object in session
  
  for (i in 1:length(x)){
    print(x[[i]])
  }
}