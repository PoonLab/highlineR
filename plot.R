summary.Data <- function(data, ...) {
  # @arg session highlineR Data object to be summarized
  # @return list of summarized variables
  
  fn <- data$path # filename
  dt <- class(data)[1] # datatype of file
  seq_num <- length(data$raw_seq) # number of sequences
  seq_len <- nchar(data$raw_seq[[1]]$sequence) # length of sequence
  var_num <- length(ls(data$compressed)) # number of variants
  
  # convert environment of variant counts to matrix
  var_counts <- (as.data.frame(as.list(data$compressed))) 
  # order variant matrix by abundance
  var_counts.sorted <- t(var_counts[order(var_counts, decreasing = TRUE)])
  colnames(var_counts.sorted) <- "freq"
  
  abun_var <- rownames(var_counts.sorted)[1] # most abundant variant
  
  res <- list(fn = fn,
              dt = dt,
              seq_num = seq_num,
              seq_len = seq_len,
              var_num = var_num,
              abun_var = abun_var,
              vars = var_counts.sorted)
  class(res) <- "summary.Data"
  res
}

print.summary.Data <- function(x, ...) {
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
  res <- eapply(session, summary)
  class(res) <- "summary.session"
  res
}

print.summary.session <- function(x, ...) {
  for (i in 1:length(x)){
    print(x[[i]])
  }
}

plot.session <- function(session = highlineR.data, master) {
  
}