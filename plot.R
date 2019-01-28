library(ggplot2)
library(ggpubr)
library(grid)

# TODO: custom arrange
# TODO: shared axes (facet? would fix ggarrange) 
# TODO: specify reading frame


plot.session <- function(session, mode, master, sort_by, ...) {
  # @arg session highlineR session of Data objects to be plotted
  
  if (missing(mode) && missing(sort_by)){
    res <- eapply(session, plot, session_plot = T)
  }
  else if(missing(mode)) {
    res <- eapply(session, plot, sort_by = sort_by, session_plot = T)
  }
  else if(missing(sort_by)) {
    res <- eapply(session, plot, mode = mode, session_plot = T)
  }
  else {
    res <- eapply(session, plot, sort_by = sort_by, mode = mode, session_plot = T)
  }
  
  figure <- ggarrange(plotlist = res, nrow = 1, ncol = length(ls(session)), common.legend = TRUE)
  annotate_figure(figure,
                  bottom = text_grob("Alignment Position", size = rel(24)))
}


plot.Data <- function(data, session_plot = F, mode = "mismatch", master = data$master, sort_by = "similarity", ...) {
  # @arg data highlineR data object to be plotted
  # @arg master sequence to be used as master, optional, default: most abundant sequence
  # @arg order ("similarity", "frequency") mode to sort sequences by
  
  if (length(data$compressed) == 0) {
    warning(paste("File", data$path, "ignored. Run highlineR::compress(...)"))
  }
  
  mode <- match.arg(mode, c("mismatch", "tvt", "svn"))
  if (inherits(data, "amino acid") && mode %in% c("tvt", "svn")) {
    warning ("Error: amino acid data cannot be run in this mode. Plotting mismatches compared to master")
    mode <- "mismatch"
  }
  
  # format data for plotting
  res <- plot_init(data, mode = mode, master = master, sort_by = sort_by)
  
  data_matrix <- res[[1]]
  rel_abun <- res[[2]]
  
  rel_abun_p <- rep(NA, length(rel_abun))
  rel_abun_p[1] <- sqrt(rel_abun[1])/2
  l <- 0
  
  for (i in 2:length(rel_abun_p)) {
    p <- sqrt(rel_abun[i-1]) + sqrt(rel_abun[i])/2 + 5
    rel_abun_p[i] <- p + l
    l <- l + p
  }

  # sequence labels for plotting
  seqs <- ls(data$compressed) # list of variants
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
  

  gg <- ggplot(data_matrix, aes(x=position, y=as.numeric(as.character(seq_plot_pos)))) +
    # plot horizontal lines using relative abundances as line height
    geom_hline(yintercept = rel_abun_p, 
              size = sqrt(as.numeric(as.character(rel_abun))), 
               color = "grey") +
    # plot vertical lines for mismatches
    geom_point(shape = "|", 
               aes(colour = value, 
                   size=sqrt(as.numeric(as.character(rel_abun))),
                   stroke = 0)) +
    scale_y_continuous(limits = c(min(rel_abun_p) - 0.1, max(rel_abun_p) + sqrt(max(rel_abun))/2), breaks=rel_abun_p, labels = seq_groups) +
    scale_x_continuous(limits = c(1, NA)) +
    labs(x = "Alignment Position", y = element_blank(), title = "Mismatches compared to master", subtitle = data$path) +
    theme_linedraw() +
    theme(axis.text.x = element_text(size = rel(2)), axis.title.x = element_text(size = rel(3))) +
    scale_size_identity()
  print(session_plot == T)
  if (session_plot == T) {
    gg <- gg + labs(x = element_blank())
  }
  if (inherits(data, "nucleotide")) {
    if (mode == "mismatch"){
      gg + scale_color_manual(name = "Legend",
                              breaks = c("A", "C", "G", "T", "-"),
                              labels = c("A", "C", "G", "T", "Gap"),
                              values = c("A" = "#00BF7D", "C" = "#00B0F6", "G" = "#A3A500", "T" = "#F8766D", "-" = "dark grey"))
    }
    else if(mode == "tvt") {
      gg + scale_color_manual(name = "Legend",
                             breaks = c("transition", "transversion"),
                             labels = c("Transition", "Transversion"),
                             values = c("transition" = "#00BF7D", "transversion" = "#00B0F6"))
    }
    else if (mode == "svn") {
      gg + scale_color_manual(name = "Legend",
                              breaks = c("synonymous", "non-synonymous"),
                              labels = c("Synonymous", "Non-Synonymous"),
                              values = c("synonymous" = "#00BF7D", "non-synonymous" = "#00B0F6"))
    }
   
  }
  else if (inherits(data, "amino acid")) {
    # TODO: custom colouring
    gg
  }
}

plot_init <- function(data, mode, master, sort_by = "similarity", ...) {
  sort_by <- match.arg(tolower(sort_by), c("similarity", "frequency"))
  
  # if master not specified, use most abundant sequence
  if (missing(master)) {
    master <- data$master
  }
  
  # identify compositional differences between sequences
  calc_seq_diff(data, mode = mode, master = master)
  
  
  
  # determine order of sequences for plotting
  seqs = NULL
  if (sort_by == "similarity"){
    seqs <- seq_simil(data$seq_diff)
  }
  else if (sort_by == "frequency"){
    seqs <- rownames(sort(data$compressed))[-1]
  }
  
  # calculate relative abundances for line thickness
  rel_abun <- calc_rel_abun(data$compressed, c(seqs, master))
  
  # format data for plotting
  data_matrix <- data_melt(data$seq_diff, seq_order = seqs, master = master)
  
  list(data_matrix, rel_abun)
}

data_melt <- function(seq_diff, seq_order, master, ...) {
  data_matrix <- data.frame(matrix(nrow = 0, ncol = 3))
  names(data_matrix) <- c("seq", "position", "value")
  for (i in 1:nrow(seq_diff)) {
    for (j in 1:ncol(seq_diff)) {
      if (!is.na(seq_diff[i, j])){
        data_matrix <- rbind(data_matrix, data.frame(seq=rownames(seq_diff)[i], position=j, value=seq_diff[i, j][[1]]))
      }
      
    }
  }
  
  # order sequences for plotting using seqs variable from above
  data_matrix$seq <- factor(data_matrix$seq, levels = c(seq_order, paste(master, "(m)")))
  
  data_matrix
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

calc_seq_diff <- function(data, mode, master) {
  # @arg data highlineR Data object to calculate sequence differences
  # @arg master sequence to which others will be compared
  
  data$seq_diff <- matrix(ncol = nchar(data$raw_seq[[1]]$sequence), nrow = length(ls(data$compressed))-1)
 
  master_seq <- strsplit(master, "")[[1]]
  
  row_num <- 1
  row_names <- NULL
  
  if (mode == "tvt") { # transitions vs transversions
    for (comp in ls(data$compressed)) { # for each sequence in environment
      if (comp != master) { # ignore master
        comp_seq <- strsplit(comp, "")[[1]]
        
        for (i in 1:length(master_seq)) { # for each positon
          if(master_seq[i] != comp_seq[i]) { # check if different
            # if A <-> G or C <-> T: transition
            if ((identical(sort(c(master_seq[i], comp_seq[i])), sort(c("A", "G")))) || (identical(sort(c(master_seq[i], comp_seq[i])), sort(c("C", "T"))))){
              data$seq_diff[row_num,i] <- "transition"
            }
            # other nucleotide substitutions: transversion
            else if (master_seq[i] != "-" && comp_seq[i] != "-") {
              data$seq_diff[row_num,i] <- "transversion"
            }
            
          }
        }
        row_names <- c(row_names, comp)
        row_num = row_num + 1
      }
    }
    rownames(data$seq_diff) <- row_names
  }
  else if(mode == "svn") { # synonymous vs non-synonymous
    for (comp in ls(data$compressed)) { # for each sequence in environment
      if (comp != master) { # ignore master
        comp_seq <- strsplit(comp, "")[[1]]
        n <- 1 # position in sequences
        
        while ((n + 2) < length(master_seq)) { # read 3 positions at a time
          master_codon <- master_seq[n:(n+2)]
          mcj <- paste(master_codon, collapse = "")
          comp_codon <- comp_seq[n:(n+2)]
          ccj <- paste(comp_codon, collapse = "")
          
          if (mcj != ccj) {
            if (identical(aaLookup(mcj), aaLookup(ccj))) {
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
    rownames(data$seq_diff) <- row_names
  }
  else if (mode == "mismatch") { # mismatches compared to master
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

aaLookup<-function(x){
  slc<-c("I","L","V","F","M","C","A","G","P","T","S","Y","W","Q","N","H","E","D","K","R","Stop")
  codon<-c("ATT, ATC, ATA","CTT, CTC, CTA, CTG, TTA, TTG","GTT, GTC, GTA, GTG","TTT, TTC","ATG","TGT, TGC",
           "GCT, GCC, GCA, GCG","GGT, GGC, GGA, GGG","CCT, CCC, CCA, CCG","ACT, ACC, ACA, ACG","TCT, TCC, TCA, TCG, AGT, AGC",
           "TAT, TAC","TGG","CAA, CAG","AAT, AAC","CAT, CAC","GAA, GAG","GAT, GAC","AAA, AAG","CGT, CGC, CGA, CGG, AGA, AGG","TAA, TAG, TGA")
  
  codon.list<-strsplit(codon,",")
  
  slc[grep(x,codon.list)]
}