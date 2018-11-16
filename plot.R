library(ggplot2)

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
              vars = var_counts.sorted)
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

sort.compressed <- function(compressed) {
  # @arg compressed environment of variant counts
  # @return matrix of variant counts sorted by abundance
  
  # convert environment of variant counts to matrix
  var_counts <- (as.data.frame(as.list(compressed))) 
  
  # order variant matrix by abundance
  var_counts.sorted <- t(var_counts[order(var_counts, decreasing = TRUE)])
  colnames(var_counts.sorted) <- "freq"
  
  var_counts.sorted
}


plot.Data <- function(data, ...) {
  data.matrix <- melt(data$seq_diff, na.rm = T)
  colnames(data.matrix) <- c("seq", "position", "value")
  ggplot(data.matrix, aes(x=position, y=seq, colour = value)) +
    geom_point(shape = "|", size=rel(10))
}

seq_diff <- function(data, master) {
  data$seq_diff <- matrix(ncol = nchar(data$raw_seq[[1]]$sequence), nrow = length(ls(data$compressed))-1)
  
  if (missing(master)) {
    master <- rownames(sort(data$compressed))[1]
  }
  master_seq <- strsplit(master, "")[[1]]
  
  row_num <- 1
  row_names <- NULL
  
  for (comp in ls(data$compressed)) { # for each sequence in environment
    if (comp != master) { # ignore master
      comp_seq <- strsplit(comp, "")[[1]]
      
      for (i in 1:length(master_seq)) { # for each positon
        if(master_seq[i] != comp_seq[i]) {
            data$seq_diff[row_num,i] <- comp_seq[i]
        }
      }
      row_names <- c(row_names, comp)
      row_num = row_num + 1
    }
  }
  rownames(data$seq_diff) <- row_names
}

plot.session <- function(session = highlineR.data, master) {
  
}
