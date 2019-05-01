#' calc_Diversity generic
#'
#' Calculates four sequence-based diversity measures in NGS Data objects.
#'
#' @param x Data object or session object containing Data objects to calculate diversity in.
#'
#' @return Returns dataframe of calculated sequence-based diversity measures.
#' @name calc_Diversity
NULL

#' @rdname calc_Diversity
#' @export
calc_Diversity <- function(x) {
  UseMethod("calc_Diversity", x)
}

#' @rdname calc_Diversity
#' @export
calc_Diversity.session <- function(x) {
  n <- ls(x)
  names(n) <- ls(x)
  res <- parallel::mclapply(n, function(data) calc_Diversity(x[[data]]), mc.cores = 20)
  dt <- data.frame(t(sapply(res,c)))
  colnames(dt) <- c("Shannon Entropy", "Percent Complexity", "Nucleotide Diversity", "Percent Diversity")
  dt
}

#' @rdname calc_Diversity
#' @export
calc_Diversity.Data <- function(x) {
  c(Shannon_Entropy(x), Percent_Complexity(x), Nucleotide_Diversity(x), Percent_Diversity(x))
}

#' @rdname calc_Diversity
#' @export
Shannon_Entropy <- function(x) {
  H <- 0 # result
  dt <- t(data.frame(strsplit(ls(x$sample), ""))) # unique variants as dataframe
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

#' @rdname calc_Diversity
#' @export
Percent_Complexity <- function(x) {
  # total number of reads
  N <- 0
  for (v in ls(x$sample)) {
    N <- N + x$sample[[v]]
  }

  # number of distinct variants divided by the total number of reads x 100
  length(x$sample)/N * 100
}

#' @rdname calc_Diversity
#' @export
Nucleotide_Diversity <- function(x) {
  pi <- 0 # result

  # number of unique variants
  n <- length(x$sample)

  # total number of reads
  N <- 0
  for (v in ls(x$sample)) {
    N <- N + x$sample[[v]]
  }

  # pairwise comparison of unique variants
  if (n > 1) {
    for (i in 2:n) {
      for (j in 1:(i-1)) {
        s1 <- strsplit(ls(x$sample)[i], "")[[1]]
        s2 <- strsplit(ls(x$sample)[j], "")[[1]]

        # ignored gapped positions
        gaps <- union(which(s1 == "-"), which(s2 == "-"))

        # length of sequence after gap removal
        l <- length(s1[-gaps])
        if (l == 0)
          next

        # number of nucleotide differences per site
        pi_ij <- (length(which(s1[-gaps] != s2[-gaps])))/l

        # sum (frequency of sequence i * frequency of sequence j * pi_ij)
        pi <- pi + ((x$sample[[paste(s1, collapse = "")]]/N) * (x$sample[[paste(s2, collapse = "")]])/N * pi_ij)
      }
    }
  }
  2 * pi
}

#' @rdname calc_Diversity
#' @export
Percent_Diversity <- function(x) {
  seqs <- vector()
  for(s in ls(x$sample)) {
    seqs <- c(seqs, rep(s, x$sample[[s]]))
  }
  y <- t(sapply(strsplit(seqs, split = ""), tolower))
  rownames(y) <-seqs

  db <- ape::as.DNAbin(y)
  res <- ape::dist.dna(db, as.matrix=T, model='raw', pairwise.deletion=T)
  mean(res[upper.tri(res)])
}
