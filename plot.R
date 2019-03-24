library(ggplot2)
library(ggpubr)

plot.session <- function(session, mode = "mismatch", master, sort_by, rf = 1, use_sample = T, ...) {
  # @arg session: highlineR session of Data objects to be plotted
  # @arg mode: mutation annotation, optional. Options: "mismatch" (default), "svn" (Synonymous versus Non-Synonymous), "tvt" (Transition versus Transversion)
  # @arg master: sequence to be used as master, optional, default: most abundant sequence
  # @arg sort_by: how should sequences be ordered in plot. Options: "similarity", "frequency".
  # @arg rf: reading frame for determining Synonymous versus Non-Synonymous mutations, optional. Options: 1 (default), 2, 3.
  # @return plot of all Data in session
  
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
  
  figure <- ggarrange(plotlist = res, common.legend = TRUE)
  annotate_figure(figure,
                  bottom = text_grob("Alignment Position", size = rel(18)),
                  top = text_grob(modetotitle(mode), size = rel(20)))
}


plot.Data <- function(data, session_plot = F, mode = "mismatch", master = data$master, sort_by = "similarity", rf = 1, use_sample = T, ...) {
  # @arg data: highlineR data object to be plotted
  # @arg session_plot: logical. is plot being produced part of a session plot. default: FALSE
  # @arg mode: mutation annotation, optional. Options: "mismatch" (default), "svn" (Synonymous versus Non-Synonymous), "tvt" (Transition versus Transversion)
  # @arg master: sequence to be used as master, optional, default: most abundant sequence
  # @arg sort_by: how should sequences be ordered in plot. Options: "similarity", "frequency".
  # @arg rf: reading frame for determining Synonymous versus Non-Synonymous mutations, optional. Options: 1 (default), 2, 3.
  # @return plot of data showing mutations in variants compared to master

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
  if (mode == "svn") {
    if (!rf %in% 1:3) {
      stop("Error: invalid reading frame selected.")
    }
    else {
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
  gg <- ggplot(data_matrix, aes(x=position, y=as.numeric(as.character(seq_plot_pos)))) +
    # plot horizontal lines using relative abundances as line height
    geom_hline(yintercept = rel_abun_p, 
              size = sqrt(as.numeric(as.character(rel_abun))), 
               color = "grey") +
    scale_x_continuous(limits = c(1, nchar(data$master))) +
    theme_classic() +
    theme(axis.text = element_text(size = rel(1))) +
    theme(legend.position = "right") +
    scale_size_identity() +
    guides(colour = guide_legend(override.aes = list(size=5)))
  
  if (nrow(data$seq_diff) > 1) {
    # plot vertical lines for mismatches
    gg <- gg +  geom_point(shape = "|", 
                           aes(colour = value, 
                               size=sqrt(as.numeric(as.character(rel_abun)))+1,
                               stroke = 0)) +
      scale_y_continuous(limits = c(min(rel_abun_p) - 0.1, max(rel_abun_p) + sqrt(max(rel_abun))/2), breaks=rel_abun_p, labels = seq_groups)
  }
  else {
    gg <- gg + scale_y_continuous(breaks=rel_abun_p, labels = seq_groups)
  }
  
  if (session_plot == T) {
    gg <- gg + labs(x = element_blank(), y = element_blank(), title = element_blank(), subtitle = filename[length(filename)])
    print("here")
  }
  else {
    gg <- gg + labs(x = "Alignment Position", y = element_blank(), title = modetotitle(mode), subtitle = filename[length(filename)]) +
      theme(axis.title.x = element_text(size = rel(2)))
  }
  
  print(mode)
  if (inherits(data, "nucleotide")) {
    if (mode == "mismatch"){
      gg <- gg + scale_color_manual(name = "Legend",
                              breaks = c("A", "C", "G", "T", "-", "del"),
                              labels = c("A", "C", "G", "T", "Gap", ""),
                              values = c("A" = "#00BF7D", "C" = "#00B0F6", "G" = "#A3A500", "T" = "#F8766D", "-" = "dark grey", "del" = "white"),
                              drop = FALSE)
    }
    else if(mode == "tvt") {
      gg <- gg + scale_color_manual(name = "Legend",
                              drop = FALSE,
                             breaks = c("transition", "transversion", "del"),
                             labels = c("Transition", "Transversion", ""),
                             values = c("transition" = "#00BF7D", "transversion" = "#00B0F6", "del" = "white"))
    }
    else if (mode == "svn") {
      gg <- gg + scale_color_manual(name = "Legend",
                              drop = FALSE,
                              breaks = c("synonymous", "non-synonymous", "del"),
                              labels = c("Synonymous", "Non-Synonymous", ""),
                              values = c("synonymous" = "#00BF7D", "non-synonymous" = "#00B0F6", "del" = "white"))
    }
   
  }
  else if (inherits(data, "amino acid")) {
  gg <- gg + scale_color_manual(name = "Legend",
                                  drop = FALSE,
                                  breaks = c("H", "DE", "KNQR", "M", "ILV", "FWY", "C", "AGST", "P", "-", "del"),
                                  labels = c("His", "Asp, Glu", "Lys, Asn, Gln, Arg", "Met", "Ile, Leu, Val", "Phe, Trp, Tyr", "Cys", "Ala, Gly, Ser, Thr", "Pro", "Gap", ""),
                                  values = c("H" = "blue", "DE" = "navyblue", "KNQR" = "skyblue", "M" = "darkgreen", "ILV" = "green", "FWY" = "magenta", "C" = "red", "AGST" = "orange", "P" = "yellow", "-" = "dark grey", "del" = "white"))
  }
  # ggsave(paste0(data$path, ".pdf"), plot = gg, device = "pdf", width = 10, height = 16, units = "in", dpi = 300)
  gg
}

plot_init <- function(data, compressed, mode, master, sort_by = "similarity", rf, ...) {
  # @arg data: highlineR data object to be plotted
  # @arg mode: mutation annotation, optional. Options: "mismatch" (default), "svn" (Synonymous versus Non-Synonymous), "tvt" (Transition versus Transversion)
  # @arg master: sequence to be used as master, optional, default: most abundant sequence
  # @arg sort_by: how should sequences be ordered in plot. Options: "similarity", "frequency".
  # @arg rf: reading frame for determining Synonymous versus Non-Synonymous mutations, optional. Options: 1 (default), 2, 3.
  # @return list of mutation data for ploting and relative abundances of variants
  
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
  seqs = NULL
  if (sort_by == "similarity"){
    seqs <- c(seq_simil(subset(data$seq_diff, rownames(data$seq_diff) != master)), master)
  }
  else if (sort_by == "frequency"){
    seqs <- rownames(sort(compressed))
  }
  
  # calculate relative abundances for line thickness
  print(".... Calculating relative frequencies")
  rel_abun <- calc_rel_abun(compressed, seqs)
  
  # format data for plotting
  print(".... Formatting data for plotting")
  data_matrix <- data_melt(data$seq_diff, seq_order = seqs, master = master)

  list(data_matrix, rel_abun)
}

calc_seq_diff <- function(data, compressed, mode, master, rf) {
  # @arg data: highlineR Data object to calculate compositional differences
  # @arg mode: mutation annotation. Options: "mismatch", "svn", "tvt"
  # @arg master sequence to which others will be compared
  # @arg rf: reading frame for determining Synonymous versus Non-Synonymous mutations, optional. Options: 1, 2, 3.
  # identifies compositional differences between sequences in data including annotations if applicable. stores result in data$seq_dif
  
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

seq_simil <- function(seq_diff) {
  # @arg seq_diff: matrix of compositional differences between sequences
  # @return vector of sequences ordered by number of mutations (similarity) compared to master
  
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

calc_rel_abun <- function(compressed, seqs) {
  # @arg compressed: environment of variant counts
  # @arg seqs: vector of sequences present in compressed in order of interest
  # @return vector of read counts of seqs from compressed in order specified
  
  res <- NULL
  for (seq in seqs) {
    # retrieve sequence read counts from compressed
    res <- c(res, get(seq, envir = compressed))
  }
  res
}

data_melt <- function(seq_diff, seq_order, master, ...) {
  # @arg seq_diff: matrix of compositional differnces between sequences
  # @arg seq_order: vector of sequences in order they should be plotted
  # @arg master: sequence to which others have been compared
  # @return matrix of compositional differences between sequences in format required for ggplot plotting
  
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

aaLookup<-function(c){
  # @arg c: sequence of three bases (codon)
  # @return single letter code of amino acid encoded by codon
  
  slc<-c("I","L","V","F","M","C","A","G","P","T","S","Y","W","Q","N","H","E","D","K","R","Stop")
  codon<-c("ATT, ATC, ATA","CTT, CTC, CTA, CTG, TTA, TTG","GTT, GTC, GTA, GTG","TTT, TTC","ATG","TGT, TGC",
           "GCT, GCC, GCA, GCG","GGT, GGC, GGA, GGG","CCT, CCC, CCA, CCG","ACT, ACC, ACA, ACG","TCT, TCC, TCA, TCG, AGT, AGC",
           "TAT, TAC","TGG","CAA, CAG","AAT, AAC","CAT, CAC","GAA, GAG","GAT, GAC","AAA, AAG","CGT, CGC, CGA, CGG, AGA, AGG","TAA, TAG, TGA")
  
  codon.list<-strsplit(codon,",")
  
  # find codon in codon.list and retrun corresponding single letter code from slc
  slc[grep(c,codon.list)]
}

modetotitle <- function(mode) {
  # @arg mode: code for mutation annotion being plotted
  # @return string representation of mode
  
  titles <- c(mismatch = "Mismatches Compared to Master",
              tvt = "Transitions and Transversions",
              svn = "Synonymous and Non-Synonymous Mutations")
  titles[[mode]]
}
