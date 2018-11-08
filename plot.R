data <- highlineR.data$`/home/lisamonique/Documents/highlineR/data/test.fa`

summary.Data <- function(data) {
  cat("Summary:\n")
  cat("Data type:", class(data)[1], "\n")
  cat("Number of sequences:", length(data$raw_seq), "\n")
  cat("Sequence length:", nchar(data$raw_seq[[1]]$sequence), "\n")
  cat("Number of variants:", length(ls(data$compressed)), "\n\n")
  
  ls.str(data$compressed)
}

plot.session <- function(session = highlineR.data, master) {
  
}