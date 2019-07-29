#' highlineR plot functions
#' 
#' Generic plot functions for highlineR objects
#' 
#' @param x Data object or session object containing Data objects to plot.
#' @param mode A character string representing the desired mutation annotation for plotting. 
#' Options: "mismatch" (default), "svn" (Synonymous versus Non-Synonymous), "tvt" (Transition 
#' versus Transversion).
#' @param master A character string representing the sequence to which other sequences 
#' should be compared. By default, the most abundant sequence is selected.
#' @param sort_by A character string representing how the sequences should be ordered in 
#' a plot. Options: "similarity", "frequency".
#' @param rf An integer specifying which reading frame should be used when determining 
#' Synonymous vs Non-Synonymous mutations. Options: 1 (default), 2, 3.
#' @param use_sample A logical value specifying which environment should be plotted. If 
#' \code{True}, then the \code{sample} environment of the Data objects is plotted. If \code{False}, the complete \code{compressed} environment is plotted.
#' 
#' @return Object of class "ggplot2"
#' @name plot
NULL

#' @rdname plot
#' @export
plot.session <- function(x, mode = "mismatch", master = NA, sort_by = NA, rf = 1, 
                         use_sample = T, quiet=F, size.title=16, size.xlab=14, ...) {
  # plot.session is a generic plot function for objects of class "session" (highlineR)
  # This function is a wrapper for calling plot.Data
  #
  # Args:
  #   x: Data object or session object containing Data objects to plot.
  #   mode: A character string representing the desired mutation annotation for 
  #         plotting.  Options: "mismatch" (default), "svn" (Synonymous versus Non-Synonymous),
  #         "tvt" (Transition versus Transversion).
  #   master: A character or integer vector representing the sequence to which other 
  #           sequences should be compared. By default, the most abundant 
  #           sequence is selected.
  #   sort_by: A character string representing how the sequences should be ordered 
  #            in a plot. Options: "similarity", "frequency".
  #   rf: An integer specifying which reading frame should be used when determining
  #       Synonymous vs Non-Synonymous mutations. Options: 1 (default), 2, 3.
  #   use_sample: A logical value specifying which environment should be plotted.
  #               If \code{True}, then the \code{sample} environment of the Data
  #               objects is plotted. If \code{False}, the complete \code{compressed}
  #               environment is plotted.
  
  # apply different master sequence arguments to each data set
  n <- length(x)
  res <- lapply(1:n, function(i) {
    plot(  # calls plot.Data
      get(ls(x)[i], x),  # extract Data object
      mode=mode,  # cannot mix modes
      master=ifelse(length(master)==n, master[i], master),
      sort_by=ifelse(length(sort_by)==n, sort_by[i], sort_by),
      rf=ifelse(length(rf)==n, rf[i], rf),
      use_sample=ifelse(length(use_sample)==n, use_sample[i], use_sample),
      session_plot=T,
      quiet=quiet
    )
  })

  #res <- eapply(x, plot, mode = mode, master = master, sort_by = sort_by, 
  #              session_plot = T, rf = rf, use_sample = use_sample, quiet=quiet)

  figure <- ggpubr::ggarrange(plotlist = rev(res), common.legend = TRUE)
  ggpubr::annotate_figure(figure,
                  bottom = ggpubr::text_grob("Alignment Position", 
                                             size = ggplot2::rel(size.xlab)),
                  top = ggpubr::text_grob(modetotitle(mode), 
                                          size = ggplot2::rel(size.title)))
}

#' @rdname plot
#' @export
plot.Data <- function(x, mode = NA, master = NA, sort_by = NA, 
                      rf = 1, use_sample = T, quiet = F, session_plot = F, ...) {
  # plot.Data converts a highlineR object of class Data into a 
  # matrix for plotting, and then runs call.ggplot to generate the 
  # plot object
  
  # default arguments
  if (is.na(mode)) mode='mismatch'
  if (is.na(master)) master = x$master
  if (is.na(sort_by)) sort_by = 'similarity'
  
  data <- x
  
  if (!quiet) print(paste("Plotting:", data$path))
  
  # check inputs
  if (use_sample == T) {
    compressed <- data$sample
  } else {
    compressed <- data$compressed
  }
  if (length(compressed) == 0) {
    stop(paste("File", data$path, "ignored. Run highlineR::compress(...)"))
  }

  mode <- match.arg(mode, c("mismatch", "tvt", "svn"))
  
  if (inherits(data, "amino acid") && mode %in% c("tvt", "svn")) {
    warning ("Error: amino acid data cannot be run in this mode. Plotting mismatches compared to master")
    mode <- "mismatch"
  }

  if (is.numeric(master)) {
    # user specified master sequence by index
    if (length(master) > 1) {
      stop("Error: only one master sequence can be selected.")
    }
    master <- data$raw_seq[[master]]$sequence
  }

  # format data for plotting
  if (!quiet) {
    print(".... Initializing Plot")
    print("........ Parameters:")
    print(paste0("............ Mode: ", mode))
    print(paste0("............ Master: ", master))
    print(paste0("............ Sort by: ", sort_by))
    print(paste0("............ Use sample?: ", use_sample))
  }

  
  if (mode == "svn") {
    # label synonymous vs. non-synonymous given reading frame <rf>
    if (!rf %in% 1:3) {
      stop("Error: invalid reading frame selected.")
    }
    else {
      if (!quiet) print(paste0("............ Reading frame: ", rf))
      res <- plot_init(data, compressed = compressed, mode = mode, master = master, 
                       sort_by = sort_by, rf = rf, quiet=quiet)
    }
  }
  else {
    res <- plot_init(data, compressed = compressed, mode = mode, master = master, 
                     sort_by = sort_by, quiet=quiet)
  }

  data_matrix <- res[[1]]
  rel_abun <- res[[2]]

  if (!quiet) print(".... Determining variant positions")
  rel_abun_p <- rep(NA, length(rel_abun))
  rel_abun_p[1] <- sqrt(rel_abun[1])/2
  l <- 0

  if (length(compressed) > 1) {
    for (i in 2:length(rel_abun_p)) {
      p <- sqrt(rel_abun[i-1]) + sqrt(rel_abun[i])/2 + 5
      rel_abun_p[i] <- p + l
      l <- l + p
    }
  }

  # sequence labels for plotting
  seqs <- ls(compressed) # list of variants
  seq_groups <- vector()
  for (s in seqs) {
    for (i in 1:length(data$raw_seq)) {
      if (data$raw_seq[[i]]$sequence == s) {
        seq_groups <- c(seq_groups, data$raw_seq[[i]]$header)
        break
      }
    }
  }
  seq_groups <- seq_groups[
    order(as.numeric(factor(
      gsub(master, paste(master, "(m)"), seqs, perl=T), 
      levels = levels(data_matrix$seq)
    )))
  ]
  seq_groups[length(seqs)] <- paste(seq_groups[length(seqs)], "(m)")

  data_matrix$seq_plot_pos <- data_matrix$seq
  levels(data_matrix$seq_plot_pos) <- rel_abun_p
  data_matrix$rel_abun <- data_matrix$seq
  levels(data_matrix$rel_abun) <- rel_abun
  
  if (mode == "mismatch") {
    if (inherits(data, "nucleotide")) {
      data_matrix$value <- factor(data_matrix$value, levels = c("A", "C", "G", "T", "-", "del"))
    }
    else if (inherits(data, "amino acid")) {
      data_matrix$value <- factor(data_matrix$value, levels = c("H", "DE", "KNQR", "M", "ILV", "FWY", "C", "AGST", "P", "-", "del"))
    }
  }
  else if (mode == "tvt") {
    data_matrix$value <- factor(data_matrix$value, levels = c("transition", "transversion", "del"))
  }
  else if (mode == "svn") {
    data_matrix$value <- factor(data_matrix$value, levels = c("synonymous", "non-synonymous", "del"))
  }

  filename <- strsplit(data$path, "/")[[1]]

  if (!quiet) print(".... Done")
  call.ggplot(data, data_matrix, seq_groups, rel_abun, rel_abun_p, filename, session_plot, mode)
}



call.ggplot <- function(data, data_matrix, seq_groups, rel_abun, rel_abun_p, filename, 
                        session_plot, mode) {
  # Use ggplot2 to generate visualization of data matrix from plot.Data
  # Internal.
  # TODO: make legend optional
  #
  # Args:
  #   data_matrix: matrix generated by plot_init
  #   seq_groups:
  #   filename: 
  #   session_plot:  
  #   mode:  how to annotate differences by colour
  # Returns:
  #   Object of class 'ggplot'
  
  gg <- ggplot2::ggplot(data_matrix, 
                        ggplot2::aes(x=position, 
                                     y=as.numeric(as.character(seq_plot_pos)))) +
    # plot horizontal lines using relative abundances as line height
    ggplot2::geom_hline(yintercept = rel_abun_p,
                        size = sqrt(as.numeric(as.character(rel_abun))),
                        color = "grey") +
    ggplot2::scale_x_continuous(limits = c(1, nchar(data$master))) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = ggplot2::rel(1))) +
    ggplot2::theme(legend.position = "right") +
    ggplot2::scale_size_identity() +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=5)))
  
  if (nrow(data$seq_diff) > 1) {
    # plot vertical lines for mismatches
    gg <- gg +  ggplot2::geom_point(shape = "|",
                                    ggplot2::aes(colour = value,
                                                 size=sqrt(as.numeric(as.character(rel_abun)))+1,
                                                 stroke = 0)) +
      ggplot2::scale_y_continuous(limits = c(min(rel_abun_p) - 0.1, max(rel_abun_p) + sqrt(max(rel_abun))/2), breaks=rel_abun_p, labels = seq_groups)
  }
  else {
    gg <- gg + ggplot2::scale_y_continuous(breaks=rel_abun_p, labels = seq_groups)
  }
  
  if (session_plot == T) {
    gg <- gg + ggplot2::labs(x = ggplot2::element_blank(), y = ggplot2::element_blank(), title = ggplot2::element_blank(), subtitle = filename[length(filename)])
  }
  else {
    gg <- gg + ggplot2::labs(x = "Alignment Position", y = ggplot2::element_blank(), title = modetotitle(mode), subtitle = filename[length(filename)]) +
      ggplot2::theme(axis.title.x = ggplot2::element_text(size = ggplot2::rel(2)))
  }
  
  if (inherits(data, "nucleotide")) {
    if (mode == "mismatch"){
      gg <- gg + ggplot2::scale_color_manual(name = "Legend",
                                             breaks = c("A", "C", "G", "T", "-", "del"),
                                             labels = c("A", "C", "G", "T", "Gap", ""),
                                             values = c("A" = "#00BF7D", "C" = "#00B0F6", "G" = "#A3A500", "T" = "#F8766D", "-" = "#646464", "del" = "white"),
                                             drop = FALSE)
    }
    else if(mode == "tvt") {
      gg <- gg + ggplot2::scale_color_manual(name = "Legend",
                                             drop = FALSE,
                                             breaks = c("transition", "transversion", "del"),
                                             labels = c("Transition", "Transversion", ""),
                                             values = c("transition" = "red", "transversion" = "blue", "del" = "white"))
    }
    else if (mode == "svn") {
      gg <- gg + ggplot2::scale_color_manual(name = "Legend",
                                             drop = FALSE,
                                             breaks = c("synonymous", "non-synonymous", "del"),
                                             labels = c("Synonymous", "Non-Synonymous", ""),
                                             values = c("synonymous" = "red", "non-synonymous" = "blue", "del" = "white"))
    }
    
  }
  else if (inherits(data, "amino acid")) {
    gg <- gg + ggplot2::scale_color_manual(name = "Legend",
                                           drop = FALSE,
                                           breaks = c("H", "DE", "KNQR", "M", "ILV", "FWY", "C", "AGST", "P", "-", "del"),
                                           labels = c("His", "Asp, Glu", "Lys, Asn, Gln, Arg", "Met", "Ile, Leu, Val", "Phe, Trp, Tyr", "Cys", "Ala, Gly, Ser, Thr", "Pro", "Gap", ""),
                                           values = c("H" = "blue", "DE" = "navyblue", "KNQR" = "skyblue", "M" = "darkgreen", "ILV" = "green", "FWY" = "magenta", "C" = "red", "AGST" = "orange", "P" = "yellow", "-" = "dark grey", "del" = "white"))
  }
  gg
}



plot_init <- function(data, compressed, mode, master, sort_by = "similarity", 
                      rf, quiet=FALSE, ...) {
  # plot_init carries out plot initialization methods
  # Internal.
  #
  # Args:
  #   data:  Object of class "Data"
  #   compressed:  Character vector of unique sequence variants
  #   mode:  One of 'mismatch', 'svn' and 'tvt'
  #   sort_by:  One of 'similarity' (default) and 'frequency'
  # Returns:
  #   A list of mutations and relative abundances of variants for plotting
  sort_by <- match.arg(tolower(sort_by), c("similarity", "frequency"))
  mode <- match.arg(tolower(mode), c("mismatch", "svn", "tvt"))

  # if master not specified, use most abundant sequence
  if (missing(master)) {
    master <- data$master
  }

  if (!quiet) print(".... Identifying mutations")
  # identify compositional differences between sequences
  calc_seq_diff(data, compressed = compressed, mode = mode, master = master, rf = rf)


  # determine order of sequences for plotting
  if (!quiet) print(".... Determining sequence order")
  seq_order <- NULL
  if (sort_by == "similarity"){
    seq_order <- c(seq_simil(subset(data$seq_diff, rownames(data$seq_diff) != master)), master)
  }
  else if (sort_by == "frequency"){
    seq_order <- rownames(sort(compressed))
    seq_order <- c(seq_order[which(seq_order != master)], master)
  }

  # calculate relative abundances for line thickness
  if (!quiet) print(".... Calculating relative frequencies")
  rel_abun <- calc_rel_abun(compressed, seq_order)

  # format data for plotting
  if (!quiet) print(".... Formatting data for plotting")
  data_matrix <- data_melt(data$seq_diff, seq_order = seq_order, master = master)

  list(data_matrix, rel_abun)
}



calc_seq_diff <- function(data, compressed, mode, master, rf) {
  # calc_seq_diff identifies the compositional differences between sequences
  # in the data, including annotations if applicable.  Results are stored in the
  # Data object's seq_diff data frame (modified in place).
  # Internal.
  #
  # Args:
  #   data:  Object of class Data (highlineR)
  #   compressed:  character vector of unique sequence variants
  #   mode:  One of 'tvt', 'svn', or 'mismatch'
  #   master:  master sequence
  #   rf:  <optional> Reading frame, required for 'svn' mode
  # Returns:
  #   An LxN matrix of compositional differences between sequences.
  #   where L columns map to sequence length and N rows map to number
  #   of sequences).  Row names contain the sequences.
  
  # initialize matrix 
  data$seq_diff <- matrix(ncol = nchar(data$raw_seq[[1]]$sequence), 
                          nrow = length(ls(compressed)))

  master_seq <- strsplit(master, "")[[1]]

  row_num <- 1
  row_names <- NULL

  # if position in master is gap, record as deletion in variant
  deletions <- which(master_seq == "-")
  data$seq_diff[nrow(data$seq_diff), deletions] <- "del"

  if (mode == "tvt") { # transitions vs transversions
    md <- which(master_seq != "-")
    
    for (comp in ls(compressed)) { # for each sequence in environment
      
      if (comp == master) { next } # ignore master sequence
    
      # convert to character vector of nucleotides
      comp_seq <- strsplit(comp, "")[[1]]

      # index the nucleotide differences
      mismatches <- which(master_seq != comp_seq)
      
      # ignore gaps in the master sequence
      valid <- intersect(md, which(comp_seq != "-"))
      mismatches <- intersect(mismatches, valid)

      for (i in mismatches) {
        # if A <-> G or C <-> T: transition
        subn <- sort(c(master_seq[i], comp_seq[i]))
        if (identical(subn, c("A", "G")) || identical(subn, c("C", "T"))) {
          data$seq_diff[row_num,i] <- "transition"
        } else {
          # potentially a transversion
          if (master_seq[i] != "-" && comp_seq[i] != "-") {
            data$seq_diff[row_num,i] <- "transversion"
          }
        }
      }

      row_names <- c(row_names, comp)
      row_num = row_num + 1
    }
    rownames(data$seq_diff) <- c(row_names, master)
  }
  
  else if(mode == "svn") { # synonymous vs non-synonymous
    for (comp in ls(compressed)) { # for each sequence in environment
      # ignore master sequence
      if (comp == master) { next }
    
      comp_seq <- strsplit(comp, "")[[1]]
      n <- rf # start position in sequences (depends on reading frame selected)

      while ((n + 2) < length(master_seq)) { # read 3 positions at a time
        master_codon <- master_seq[n:(n+2)] # codon as character vector
        mcj <- paste(master_codon, collapse = "") # codon as string
        comp_codon <- comp_seq[n:(n+2)]
        ccj <- paste(comp_codon, collapse = "")

        if (mcj != ccj) { # check if codons different
          if ("-" %in% master_codon || "-" %in% comp_codon) {
            n = n + 3
            next
          }
          if (identical(aaLookup(mcj), aaLookup(ccj))) { # check if encoded amino acid different
            data$seq_diff[row_num, (n - 1 + which(master_codon != comp_codon))] <- "synonymous"
          }
          else{
            data$seq_diff[row_num,(n - 1 + which(master_codon != comp_codon))] <- "non-synonymous"
          }
        }
        n = n + 3
      }
      row_names <- c(row_names, comp)
      row_num = row_num + 1
    }
    rownames(data$seq_diff) <- c(row_names, master)
  }
  
  else if (mode == "mismatch") { # mismatches compared to master
    
    if (inherits(data, "nucleotide")) {
      for (comp in ls(compressed)) { # for each sequence in environment
        if (comp == master) { next }  # ignore master
        comp_seq <- strsplit(comp, "")[[1]]

        # record mismatches
        mismatches <- which(master_seq != comp_seq)
        data$seq_diff[row_num, mismatches] <- comp_seq[mismatches]

        # if position in master is gap, do not record mutation in variant
        # data$seq_diff[row_num, intersect(mismatches, deletions)] <- NA

        row_names <- c(row_names, comp)
        row_num = row_num + 1
      }
      rownames(data$seq_diff) <- c(row_names, master)
    }
    
    else if (inherits(data, "amino acid")) {
      for (comp in ls(compressed)) { # for each sequence in environment
        if (comp == master) { next } # ignore master
        comp_seq <- strsplit(comp, "")[[1]]

        # record mismatches
        mismatches <- which(master_seq != comp_seq)
        data$seq_diff[row_num, mismatches] <- comp_seq[mismatches]

        # if position in master is gap, do not record mutation in variant
        data$seq_diff[row_num, intersect(mismatches, deletions)] <- NA

        row_names <- c(row_names, comp)
        row_num = row_num + 1
      }
      
      rownames(data$seq_diff) <- c(row_names, master)
      data$seq_diff[which(data$seq_diff %in% c("D", "E"))] <- "DE"
      data$seq_diff[which(data$seq_diff %in% c("I", "L", "V"))] <- "ILV"
      data$seq_diff[which(data$seq_diff %in% c("F", "W", "Y"))] <- "FWY"
      data$seq_diff[which(data$seq_diff %in% c("A", "G", "S", "T"))] <- "AGST"
      data$seq_diff[which(data$seq_diff %in% c("K", "N", "Q", "R"))] <- "KNQR"
    }

  }
}




seq_simil <- function(seq_diff) {
  # seq_simil orders sequences based on similarity
  # Internal.
  # 
  # Args:
  #   seq_diff:  Return value of calc_seq_diff (character matrix)
  # Returns:
  #   A character vector of sequences in ascending order with respect 
  #   to the number of mutations away from master sequence.
  if (length(rownames(seq_diff)) > 1) {
    simil <- NULL
    
    # iterate over sequences (contained as row names)
    for (r in rownames(seq_diff)) {
      # non-missing entries in row indicate differences from master seq
      simil <- c(simil, length(which(!is.na(seq_diff[r,]))))
    }

    # return sequences ordered by similarity
    rownames(seq_diff)[rev(order(simil))]
  }
  else {
    rownames(seq_diff)
  }
}


sort.compressed <- function(x, decreasing, ...) {
  # Generic sort function for "compressed" environments
  # Orders sequences based on relative abundance.
  # 
  # Args:
  #   x: compressed environment object
  #   decreasing: defaults to FALSE (ascending order)
  # Returns:
  #   A Nx1 integer data frame that rearranges the provided sequences 
  #   with respect to relative abundance, where N is the number of 
  #   sequence variants.  Rownames are the sequences.
  
  # convert environment of variant counts to matrix
  var_counts <- unlist(as.list(x))

  # order variant matrix by abundance
  var_counts.sorted <- data.frame(var_counts[order(var_counts, decreasing = F)])
  colnames(var_counts.sorted) <- "freq"

  var_counts.sorted
}


calc_rel_abun <- function(compressed, seq_order) {
  # calc_rel_abun calculates sequence relative abundances.
  #
  # Args:
  #   compressed: environment containing unique sequence variants
  #   seq_order: integer vector to specify ordering of sequences in returned vector
  # Returns:
  #   Vector of read counts 
  res <- NULL
  for (seq in seq_order) {
    # retrieve sequence read counts from compressed
    res <- c(res, get(seq, envir = compressed))
  }
  res
}


data_melt <- function(seq_diff, seq_order, master, ...) {
  # Melts data matrix of sequence differences for easy plotting.
  # Internal.
  # 
  # Args:
  #   seq_diff: Objected returned from calc_seq_diff(), one row per sequence
  # Returns:
  #   Matrix of compositional differences between sequences
  
  # Initialize results matrix
  data_matrix <- data.frame(matrix(nrow = 0, ncol = 3))
  names(data_matrix) <- c("seq", "position", "value")
  
  if (nrow(seq_diff) == 1) {
    # only a single variant
    data_matrix[1,] <- c(rownames(seq_diff)[1], ncol(seq_diff), NA)
  }
  else {
    for (i in 1:nrow(seq_diff)) {
      # iterate over positions
      for (j in 1:ncol(seq_diff)) {
        if (!is.na(seq_diff[i, j])) {
          # convert data.frame to 3 columns: seq, position, value
          data_matrix <- rbind(
            data_matrix, 
            data.frame(
              seq=rownames(seq_diff)[i], 
              position=j, 
              value=seq_diff[i, j][[1]]
            )
          )
        }
      }
    }
  }

  # order sequences for plotting using seqs variable from above
  data_matrix$seq <- factor(data_matrix$seq, levels = seq_order)
  levels(data_matrix$seq)[levels(data_matrix$seq) == master] <- paste(master, "(m)")

  data_matrix
}



aaLookup <- function(codon){
  # aaLookup decodes a codon triplet into a single-letter amino acid symbol
  # TODO: does not handle mixtures
  #
  # Args:
  #   codon:  triplet to lookup
  # Returns:
  #   amino acid
  slc <- c(
    "I", "L", "V", "F", "M", "C", "A", 
    "G", "P", "T", "S", "Y", "W", "Q", 
    "N", "H", "E", "D", "K", "R", "Stop"
  )
  codons <- c(
    "ATT, ATC, ATA",  # I
    "CTT, CTC, CTA, CTG, TTA, TTG",  # L
    "GTT, GTC, GTA, GTG",  # V
    "TTT, TTC",  # F
    "ATG",  # M
    "TGT, TGC",  # C
    "GCT, GCC, GCA, GCG",  # A
    
    "GGT, GGC, GGA, GGG",  # G
    "CCT, CCC, CCA, CCG",  # P
    "ACT, ACC, ACA, ACG",  # T
    "TCT, TCC, TCA, TCG, AGT, AGC",  # S
    "TAT, TAC",  # Y
    "TGG",  # W
    "CAA, CAG",  # Q
    
    "AAT, AAC",  # N
    "CAT, CAC",  # H
    "GAA, GAG",  # E
    "GAT, GAC",  # D
    "AAA, AAG",  # K
    "CGT, CGC, CGA, CGG, AGA, AGG",  # R
    "TAA, TAG, TGA"  # Stop
  )

  codon.list<-strsplit(codons, ", ")

  # find codon in codon.list and retrun corresponding single letter code from slc
  slc[grep(codon, codon.list)]
}


# \code{modetotitle} converts internal mode code into verbose format for plotting.
# \code{modetotitle} returns string description of the mode.
modetotitle <- function(mode) {
  titles <- c(mismatch = "Mismatches Compared to Master",
              tvt = "Transitions and Transversions",
              svn = "Synonymous and Non-Synonymous Mutations")
  titles[[mode]]
}

utils::globalVariables(c("position", "seq_plot_pos", "value", "plot"))
