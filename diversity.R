library(ape)
library(parallel)

calc_Diversity <- function(x) {
  UseMethod("calc_Diversity", x)
}
calc_Diversity.session <- function(session) {
  n <- ls(session)
  names(n) <- ls(session)
  res <- mclapply(n, function(data) calc_Diversity(session[[data]]), mc.cores = 20)
  dt <- data.frame(t(sapply(res,c)))
  colnames(dt) <- c("Shannon Entropy", "Percent Complexity", "Nucleotide Diversity", "Percent Diversity")
  write.csv(dt, file = "diversity.csv")
  dt
}

calc_Diversity.Data <- function(data) {
  c(Shannon_Entropy(data), Percent_Complexity(data), Nucleotide_Diversity(data), Percent_Diversity(data))
}

Shannon_Entropy <- function(data) {
  H <- 0 # result
  dt <- t(data.frame(strsplit(ls(data$sample), ""))) # unique variants as dataframe
  rownames(dt) <- NULL
  L <- ncol(dt) # length of sequence
  n <- nrow(dt) # number of sequences
  
  # calculate entropy at each position
  for (l in 1:L) {
    v <- as.data.frame(table(dt[, l]))[, 2]/n # relative freq of each nucleotide at position l
    H <- H + -sum(v * log(v))
  }
  H/L
 }

Percent_Complexity <- function(data) {
  # total number of reads
  N <- 0
  for (v in ls(data$sample)) {
    N <- N + data$sample[[v]]
  }
  
  # number of distinct variants divided by the total number of reads x 100
  length(data$sample)/N * 100
}

Nucleotide_Diversity <- function(data) {
  pi <- 0 # result
  
  # number of unique variants
  n <- length(data$sample) 
  
  # total number of reads
  N <- 0 
  for (v in ls(data$sample)) {
    N <- N + data$sample[[v]]
  }
  
  # pairwise comparison of unique variants
  if (n > 1) {
    for (i in 2:n) {
      for (j in 1:(i-1)) {
        s1 <- strsplit(ls(data$sample)[i], "")[[1]]
        s2 <- strsplit(ls(data$sample)[j], "")[[1]]
        
        # ignored gapped positions
        gaps <- union(which(s1 == "-"), which(s2 == "-"))
        
        # length of sequence after gap removal
        l <- length(s1[-gaps]) 
        if (l == 0)
          next 
        
        # number of nucleotide differences per site
        pi_ij <- (length(which(s1[-gaps] != s2[-gaps])))/l 
        
        # sum (frequency of sequence i * frequency of sequence j * pi_ij)
        pi <- pi + ((data$sample[[paste(s1, collapse = "")]]/N) * (data$sample[[paste(s2, collapse = "")]])/N * pi_ij)
      }
    }
  }
  2 * pi
}

Percent_Diversity <- function(data) {
  # TODO: average pairwise genetic distance between two sequences, others used MEGA
  seqs <- vector()
  for(s in ls(data$sample)) {
    seqs <- c(seqs, rep(s, data$sample[[s]]))
  }
  y <- t(sapply(strsplit(seqs, split = ""), tolower))
  rownames(y) <-seqs
  
  db <- as.DNAbin(y)
  res <- dist.dna(db, as.matrix=T, model='raw', pairwise.deletion=T)
  mean(res[upper.tri(res)])
}
