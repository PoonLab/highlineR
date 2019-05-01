#' @export
summary.Data <- function(data, ...) {
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

#' @export
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

#' @export
summary.session <- function(session, ...) {
  res <- eapply(session, summary)
  class(res) <- "summary.session"
  res
}

#' @export
print.summary.session <- function(x, ...) {
  # print summary of each Data object in session

  for (i in 1:length(x)){
    print(x[[i]])
  }
}
