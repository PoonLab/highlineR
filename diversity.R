calc_Diversity <- function(x) {
  UseMethod("calc_Diversity", x)
}
calc_Diversity.session <- function(session) {
  # eapply(session, calc_Diversity)
  dt <- data.frame(matrix(ncol = 3))
  colnames(dt) <- c("Shannon Entropy", "Percent Complexity", "Nucleotide Diversity")
  for (data in ls(session)) {
    dt <- cbind(dt, calc_Diversity(session[[data]]))
  }
  dt
}

calc_Diversity.Data <- function(data) {
  print(c(Shannon_Entropy(data), Percent_Complexity(data), Nucleotide_Diversity(data)))
}

Shannon_Entropy <- function(data) {
  H <- 0 # result
  dt <- t(data.frame(strsplit(ls(data$compressed), ""))) # unique variants as dataframe
  rownames(dt) <- NULL
  L <- ncol(dt) # length of sequence
  
  # calculate entropy at each position
  for (l in 1:L) {
    v <- as.data.frame(table(dt[, l]))[, 2] # freq of each nucleotide at position l
    H <- H + -sum(v * log2(v))
  }
  H/L
 }

Percent_Complexity <- function(data) {
  # total number of reads
  N <- 0
  for (v in ls(data$compressed)) {
    N <- N + data$compressed[[v]]
  }
  
  # number of distinct variants divided by the total number of reads x 100
  length(data$compressed)/N * 100
}

Nucleotide_Diversity <- function(data) {
  pi <- 0 # result
  
  # number of unique variants
  n <- length(data$compressed) 
  
  # total number of reads
  N <- 0 
  for (v in ls(data$compressed)) {
    N <- N + data$compressed[[v]]
  }
  
  # pairwise comparison of unique variants
  for (i in 2:n) {
    for (j in 1:(i-1)) {
      s1 <- strsplit(ls(data$compressed)[i], "")[[1]]
      s2 <- strsplit(ls(data$compressed)[j], "")[[1]]
      
      # ignored gapped positions
      gaps <- union(which(s1 == "-"), which(s2 == "-"))
      
      # length of sequence after gap removal
      l <- length(s1[-gaps]) 
      if (l == 0)
        next 
      
      # number of nucleotide differences per site
      pi_ij <- (length(which(s1[-gaps] != s2[-gaps])))/l 
      
      # sum (frequency of sequence i * frequency of sequence j * pi_ij)
      pi <- pi + ((data$compressed[[paste(s1, collapse = "")]]/N) * (data$compressed[[paste(s2, collapse = "")]])/N * pi_ij)
    }
  }
  2 * pi
}

Percent_Diversity <- function(data) {
  # TODO: average pairwise genetic distance between two sequences, others used MEGA
}
