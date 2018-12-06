library(ggplot2)
library(reshape2)
library(ggpubr)

# TODO: custom melt
# TODO: custom arrange
# TODO: shared axes (facet? would fix ggarrange) 
# TODO: aa -> mutation classes
# TODO: sort by freq vs similarity

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


plot.Data <- function(data, master, order = "similarity",...) {
  # @arg data highlineR data object to be plotted
  # @arg master sequence to be used as master, optional, default: most abundant sequence
  # @arg order ("similarity", "frequency") method to sort sequences by
  
  order <- match.arg(tolower(order), c("similarity", "frequency"))
  
  # if master not specified, identify most abundant sequence
  if (missing(master)) {
    master <- rownames(sort(data$compressed))[1]
  }
  
  # identify mismatches from master
  calc_seq_diff(data, master)
  
  # order of sequences
  seqs = NULL
  if (order == "similarity"){
    seqs <- seq_simil(data$seq_diff)
  }
  else if (order == "frequency"){
    seqs <- rownames(sort(data$compressed))[-1]
  }
  
  # calculate relative abundances for line thickness
  rel_abun <- calc_rel_abun(data$compressed, c(seqs, master))
  
  # melt data for ggplot
  data_matrix <- melt(data$seq_diff, na.rm = T)
  colnames(data_matrix) <- c("seq", "position", "value")
  
  # order sequences for plotting using seqs variable from above
  data_matrix$seq <- factor(data_matrix$seq, levels = c(seqs, paste(master, "(m)")))
  
  gg <- ggplot(data_matrix, aes(x=position, y=seq, colour = value)) +
    # plot horizontal lines with relative abundances
    geom_hline(yintercept = 1:length(rel_abun), size = rel(rel_abun), color = "grey") +
    # plot vertical lines for mismatches
    geom_point(shape = "|", size=5) +
    # prevent dropping of master sequence despite no points
    scale_y_discrete(drop = FALSE) +
    labs(x = "Alignment Position", y = element_blank(), title = "Mismatches compared to master")
  gg
}

calc_seq_diff <- function(data, master) {
  # @arg data highlineR Data object to calculate sequence differences
  # @arg master sequence to which others will be compared
  
  data$seq_diff <- matrix(ncol = nchar(data$raw_seq[[1]]$sequence), nrow = length(ls(data$compressed))-1)
  
  master_seq <- strsplit(master, "")[[1]]
  
  row_num <- 1
  row_names <- NULL
  
  for (comp in ls(data$compressed)) { # for each sequence in environment
    if (comp != master) { # ignore master
      comp_seq <- strsplit(comp, "")[[1]]
      
      for (i in 1:length(master_seq)) { # for each positon
        if(master_seq[i] != comp_seq[i]) { # check if different
            data$seq_diff[row_num,i] <- comp_seq[i]
        }
      }
      row_names <- c(row_names, comp)
      row_num = row_num + 1
    }
  }
  rownames(data$seq_diff) <- row_names
  
}

calc_rel_abun <- function(compressed, seqs) {
  # @arg compressed environment of variant counts
  # @arg seqs sequence strings
  
  res <- NULL
  for (seq in seqs) {
    res <- c(res, get(seq, envir = compressed)+1)
  }
  res
}

seq_simil <- function(seq_diff) {
  # @arg seq_diff matrix of sequence differences
  
  simil <- NULL
  for (r in rownames(seq_diff)) {
    simil <- c(simil, length(which(!is.na(seq_diff[r,]))))
  }
  
  # return sequences ordered by similarity
  rownames(seq_diff)[rev(order(simil))]
}

plot.session <- function(session, master, ...) {
  # @arg session highlineR session of Data objects to be plotted
  
  res <- eapply(session, plot)
  ggarrange(plotlist = res, nrow = 1, ncol = length(ls(session)), common.legend = TRUE)
}
