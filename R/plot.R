#' Plot generic
#'
#' Plots highlineR Data objects
#'
#' @param x Data object or session object containing Data objects to plot.
#' @param mode A character string representing the desired mutation annotation for plotting. Options: "mismatch" (default), "svn" (Synonymous versus Non-Synonymous), "tvt" (Transition versus Transversion).
#' @param master A character string representing the sequence to which other sequences should be compared. By default, the most abundant sequence is selected.
#' @param sort_by A character string representing how the sequences should be ordered in a plot. Options: "similarity", "frequency".
#' @param rf An integer specifying which reading frame should be used when determining Synonymous vs Non-Synonymous mutations. Options: 1 (default), 2, 3.
#' @param use_sample A logical value specifying which environment should be plotted. If \code{True}, then the \code{sample} environment of the Data objects is plotted. If \code{False}, the complete \code{compressed} environment is plotted.
#' @param compressed Compressed environment of sequences for plotting.
#' @param seq_diff matrix of compositional differences between sequences.
#' @param seq_order string vector of sequences present in compressed environment ordered as desired for plotting.
#' @param c three-character string representing a nucleotide codon for decoding.
#'
#' @return Plots NGS data
#' @name plot
NULL

#' @rdname plot
#' @export
plot.session <- function(session, mode = "mismatch", master, sort_by, rf = 1, use_sample = T, ...) {
  if (missing(master) && missing(sort_by)){
    res <- eapply(session, plot, mode = mode, session_plot = T, rf = rf, use_sample = use_sample)
  }
  else if(missing(master)) {
    res <- eapply(session, plot, mode = mode, sort_by = sort_by, session_plot = T, rf = rf, use_sample = use_sample)
  }
  else if(missing(sort_by)) {
    res <- eapply(session, plot, mode = mode, master = master, session_plot = T, rf = rf, use_sample = use_sample)
  }
  else {
    res <- eapply(session, plot, mode = mode, master = master, sort_by = sort_by, session_plot = T, rf = rf, use_sample = use_sample)
  }

  figure <- ggpubr::ggarrange(plotlist = rev(res), common.legend = TRUE)
  ggpubr::annotate_figure(figure,
                  bottom = ggpubr::text_grob("Alignment Position", size = ggplot2::rel(18)),
                  top = ggpubr::text_grob(modetotitle(mode), size = ggplot2::rel(20)))
}

#' @rdname plot
#' @export
plot.Data <- function(data, session_plot = F, mode = "mismatch", master = data$master, sort_by = "similarity", rf = 1, use_sample = T, ...) {
  print(paste("Plotting:", data$path))
  if (use_sample == T) {
    compressed <- data$sample
  }
  else {
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
    if (length(master) > 1) {
      stop("Error: only one master sequence can be selected.")
    }
    master <- data$raw_seq[[master]]$sequence
  }

  # format data for plotting
  print(".... Initializing Plot")
  print("........ Parameters:")
  print(paste0("............ Mode: ", mode))
  print(paste0("............ Master: ", master))
  print(paste0("............ Sort by: ", sort_by))
  print(paste0("............ Use sample?: ", use_sample))
  if (mode == "svn") {
    if (!rf %in% 1:3) {
      stop("Error: invalid reading frame selected.")
    }
    else {
      print(paste0("............ Reading frame: ", rf))
      res <- plot_init(data, compressed = compressed, mode = mode, master = master, sort_by = sort_by, rf = rf)
    }
  }
  else {
    res <- plot_init(data, compressed = compressed, mode = mode, master = master, sort_by = sort_by)
  }

  data_matrix <- res[[1]]
  rel_abun <- res[[2]]

  print(".... Determining variant positions")
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
  seq_groups <- seq_groups[order(as.numeric(factor(gsub(master, paste(master, "(m)"), seqs), levels = levels(data_matrix$seq))))]
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

  print(".... Done")
  gg <- ggplot2::ggplot(data_matrix, ggplot2::aes(x=position, y=as.numeric(as.character(seq_plot_pos)))) +
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
  ggplot2::ggsave(paste0(data$path,"_", mode, ".pdf"), plot = gg, device = "pdf", width = 10, height = 12, units = "in", dpi = 300)
  gg
}

#' @rdname plot
#' @description \code{plot_init} performs plot initialization
#' @return \code{plot_init} returns a list of mutations and relative abundances of variants for plotting.
#' @keywords internal
plot_init <- function(data, compressed, mode, master, sort_by = "similarity", rf, ...) {
  sort_by <- match.arg(tolower(sort_by), c("similarity", "frequency"))
  mode <- match.arg(tolower(mode), c("mismatch", "svn", "tvt"))

  # if master not specified, use most abundant sequence
  if (missing(master)) {
    master <- data$master
  }

  print(".... Identifying mutations")
  # identify compositional differences between sequences
  calc_seq_diff(data, compressed = compressed, mode = mode, master = master, rf = rf)


  # determine order of sequences for plotting
  print(".... Determining sequence order")
  seq_order <- NULL
  if (sort_by == "similarity"){
    seq_order <- c(seq_simil(subset(data$seq_diff, rownames(data$seq_diff) != master)), master)
  }
  else if (sort_by == "frequency"){
    seq_order <- rownames(sort(compressed))
  }

  # calculate relative abundances for line thickness
  print(".... Calculating relative frequencies")
  rel_abun <- calc_rel_abun(compressed, seq_order)

  # format data for plotting
  print(".... Formatting data for plotting")
  data_matrix <- data_melt(data$seq_diff, seq_order = seq_order, master = master)

  list(data_matrix, rel_abun)
}

#' @rdname plot
#' @description \code{calc_seq_diff} identifies sequence differences.
#' @return \code{calc_seq_diff} identifies compositional differences between sequences in data including annotations if applicable and stores the result in the Data object's \code{seq_diff} dataframe.
#' @keywords internal
calc_seq_diff <- function(data, compressed, mode, master, rf) {
  data$seq_diff <- matrix(ncol = nchar(data$raw_seq[[1]]$sequence), nrow = length(ls(compressed)))

  master_seq <- strsplit(master, "")[[1]]

  row_num <- 1
  row_names <- NULL

  # if position in master is gap, record as deletion in variant
  deletions <- which(master_seq == "-")
  data$seq_diff[nrow(data$seq_diff), deletions] <- "del"

  if (mode == "tvt") { # transitions vs transversions
    md <- which(master_seq != "-")
    for (comp in ls(compressed)) { # for each sequence in environment
      if (comp != master) { # ignore master
        comp_seq <- strsplit(comp, "")[[1]]

        mismatches <- which(master_seq != comp_seq)
        valid <- intersect(md, which(comp_seq != "-"))
        mismatches <- intersect(mismatches, valid)

        for (i in mismatches) {
          # if A <-> G or C <-> T: transition
          if ((identical(sort(c(master_seq[i], comp_seq[i])), sort(c("A", "G")))) || (identical(sort(c(master_seq[i], comp_seq[i])), sort(c("C", "T"))))){
            data$seq_diff[row_num,i] <- "transition"
          }
          # other nucleotide substitutions: transversion
          else if (master_seq[i] != "-" && comp_seq[i] != "-") {
            data$seq_diff[row_num,i] <- "transversion"
          }
        }

        row_names <- c(row_names, comp)
        row_num = row_num + 1
      }
    }
    rownames(data$seq_diff) <- c(row_names, master)
  }
  else if(mode == "svn") { # synonymous vs non-synonymous
    for (comp in ls(compressed)) { # for each sequence in environment
      if (comp != master) { # ignore master
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
    }
    rownames(data$seq_diff) <- c(row_names, master)
  }
  else if (mode == "mismatch") { # mismatches compared to master
    if (inherits(data, "nucleotide")) {
      for (comp in ls(compressed)) { # for each sequence in environment
        if (comp != master) { # ignore master
          comp_seq <- strsplit(comp, "")[[1]]

          # record mismatches
          mismatches <- which(master_seq != comp_seq)
          data$seq_diff[row_num, mismatches] <- comp_seq[mismatches]

          # if position in master is gap, do not record mutation in variant
          # data$seq_diff[row_num, intersect(mismatches, deletions)] <- NA

          row_names <- c(row_names, comp)
          row_num = row_num + 1
        }
      }
      rownames(data$seq_diff) <- c(row_names, master)
    }
    else if (inherits(data, "amino acid")) {
      for (comp in ls(compressed)) { # for each sequence in environment
        if (comp != master) { # ignore master
          comp_seq <- strsplit(comp, "")[[1]]

          # record mismatches
          mismatches <- which(master_seq != comp_seq)
          data$seq_diff[row_num, mismatches] <- comp_seq[mismatches]

          # if position in master is gap, do not record mutation in variant
          data$seq_diff[row_num, intersect(mismatches, deletions)] <- NA

          row_names <- c(row_names, comp)
          row_num = row_num + 1
        }
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

#' @rdname plot
#' @description \code{seq_simil} orders sequences based on similarity.
#' @return \code{seq_simil} returns a permutation (an integer vector) which rearranges the provided sequences in ascending order of number of mutations
#' @keywords internal
seq_simil <- function(seq_diff) {
  if (length(rownames(seq_diff)) > 1) {
    simil <- NULL
    for (r in rownames(seq_diff)) {
      simil <- c(simil, length(which(!is.na(seq_diff[r,]))))
    }

    # return sequences ordered by similarity
    rownames(seq_diff)[rev(order(simil))]
  }
  else {
    rownames(seq_diff)
  }
}

#' @rdname plot
#' @description \code{sort} orders sequences based on relative abundance
#' @return \code{sort} returns a (an integer dataframe) which rearranges the provided sequences in ascending order of relative abundance.
#' @keywords internal
#' @export
sort.compressed <- function(compressed) {
  # @arg compressed environment of variant counts
  # @return matrix of variant counts sorted by abundance

  # convert environment of variant counts to matrix
  var_counts <- unlist(as.list(compressed))

  # order variant matrix by abundance
  var_counts.sorted <- data.frame(var_counts[order(var_counts, decreasing = F)])
  colnames(var_counts.sorted) <- "freq"

  var_counts.sorted
}

#' @rdname plot
#' @description \code{calc_rel_abun} calculates sequence relative abundances.
#' @return \code{calc_rel_abun} returns a vector of read counts of the sequences from the \code{compressed} environment in the order specified by \code{seqs_order}.
#' @keywords internal
calc_rel_abun <- function(compressed, seq_order) {
  res <- NULL
  for (seq in seq_order) {
    # retrieve sequence read counts from compressed
    res <- c(res, get(seq, envir = compressed))
  }
  res
}

#' @rdname plot
#' @description \code{data_melt} melts data for easy plotting.
#' @return \code{data_melt} returns a matrix of compositional differences between sequences in a format suitable for ggplot plotting.
#' @keywords internal
data_melt <- function(seq_diff, seq_order, master, ...) {
  data_matrix <- data.frame(matrix(nrow = 0, ncol = 3))
  names(data_matrix) <- c("seq", "position", "value")
  if (nrow(seq_diff) == 1) {
    data_matrix[1,] <- c(rownames(seq_diff)[1], ncol(seq_diff), NA)
  }
  else{
    for (i in 1:nrow(seq_diff)) {
      for (j in 1:ncol(seq_diff)) {
        if (!is.na(seq_diff[i, j])){
          # convert data.frame to 3 columns: seq, position, value
          data_matrix <- rbind(data_matrix, data.frame(seq=rownames(seq_diff)[i], position=j, value=seq_diff[i, j][[1]]))
        }
      }
    }
  }

  # order sequences for plotting using seqs variable from above
  data_matrix$seq <- factor(data_matrix$seq, levels = seq_order)
  levels(data_matrix$seq)[levels(data_matrix$seq) == master] <- paste(master, "(m)")

  data_matrix
}

#' @rdname plot
#' @description \code{aaLookup} decodes triplet codon into single letter amino acid.
#' @return \code{aaLookup} returns the single character code of the amino acid encoded by the codon \code{c}.
#' @keywords internal
aaLookup<-function(c){
  slc<-c("I","L","V","F","M","C","A","G","P","T","S","Y","W","Q","N","H","E","D","K","R","Stop")
  codon<-c("ATT, ATC, ATA","CTT, CTC, CTA, CTG, TTA, TTG","GTT, GTC, GTA, GTG","TTT, TTC","ATG","TGT, TGC",
           "GCT, GCC, GCA, GCG","GGT, GGC, GGA, GGG","CCT, CCC, CCA, CCG","ACT, ACC, ACA, ACG","TCT, TCC, TCA, TCG, AGT, AGC",
           "TAT, TAC","TGG","CAA, CAG","AAT, AAC","CAT, CAC","GAA, GAG","GAT, GAC","AAA, AAG","CGT, CGC, CGA, CGG, AGA, AGG","TAA, TAG, TGA")

  codon.list<-strsplit(codon,",")

  # find codon in codon.list and retrun corresponding single letter code from slc
  slc[grep(c,codon.list)]
}

#' @rdname plot
#' @description \code{modetotitle} converts internal mode code into verbose format for plotting.
#' @return \code{modetotitle} returns string description of the mode.
#' @keywords internal
modetotitle <- function(mode) {
  titles <- c(mismatch = "Mismatches Compared to Master",
              tvt = "Transitions and Transversions",
              svn = "Synonymous and Non-Synonymous Mutations")
  titles[[mode]]
}
