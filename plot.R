library(ggplot2)
library(ggpubr)
library(grid)

# TODO: custom arrange
# TODO: shared axes (facet? would fix ggarrange) 
# TODO: aa -> mutation classes
# TODO: store sequence groups
# TODO: single variant files

plot.session <- function(session, master, sort_by = "similarity", ...) {
  # @arg session highlineR session of Data objects to be plotted
  
  # facets: issues, dropping unused without dropping master, using re_abun
  res <- eapply(session, plot)
  res
  
  # data_matrix <- data.frame(matrix(nrow = 0, ncol = 4))
  # names(data_matrix) <- c("seq", "position", "value", "file")
  # 
  # for (i in 1:length(res)){
  #   file_data <- res[[i]][[1]]
  #   file_data$file <- names(res)[i]
  #   data_matrix <- rbind(data_matrix, file_data)
  # }
  # 
  # gg <- ggplot(data_matrix, aes(x=position, y=seq, colour = value)) +
  #   facet_wrap(~ file, scales = "free") +
  #   # plot horizontal lines with relative abundances
  #   # geom_hline(yintercept = 1:length(rel_abun), size = rel(rel_abun), color = "grey") +
  #   # plot vertical lines for mismatches
  #   geom_point(shape = "|", size=5) +
  #   # prevent dropping of master sequence despite no points
  #   # scale_y_discrete(drop = FALSE) +
  #   labs(x = "Alignment Position", y = element_blank(), title = "Mismatches compared to master")
  # gg
  # 
  # # grid doesn't solve anything
  # res <- eapply(session, plot)
  # grobs <- lapply(res, ggplotGrob)
  # plots <- grobs[[1]]
  # for (i in 2:length(grobs)){
  #   plots <- rbind(plots, grobs[[i]], size = "last")
  # }
  # grid.newpage()
  # grid.draw(plots)
  # 
  # ggarrange(plotlist = res, nrow = 1, ncol = length(ls(session)), common.legend = TRUE)
}

plot_init <- function(data, master, sort_by = "similarity", ...) {
  sort_by <- match.arg(tolower(sort_by), c("similarity", "frequency"))
  
  # if master not specified, use most abundant sequence
  if (missing(master)) {
    master <- data$master
  }
  
  # identify mismatches from master
  calc_seq_diff(data, master)
  
  # order of sequences
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

plot.Data <- function(data, master = data$master, sort_by = "similarity", ...) {
  # @arg data highlineR data object to be plotted
  # @arg master sequence to be used as master, optional, default: most abundant sequence
  # @arg order ("similarity", "frequency") method to sort sequences by
  
  if (length(data$compressed) == 0) {
    warning(paste("File", data$path, "ignored. Run highlineR::compress(...)"))
  }
  
  # format data for plotting
  res <- plot_init(data, master, sort_by)
  data_matrix <- res[[1]]
  rel_abun <- res[[2]]
  
  rel_abun_sum <- sqrt(rel_abun)
  for (i in 2:length(rel_abun_sum)) {
    rel_abun_sum[i] <- rel_abun_sum[i] + rel_abun_sum[i-1]
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
  print(seq_groups)
  # seq_groups <- paste0("v", 0:(length(seqs)-1)) # list of labels of form "v_n"
  seq_groups <- seq_groups[order(as.numeric(factor(gsub(master, paste(master, "(m)"), seqs), levels = levels(data_matrix$seq))))]
  seq_groups[length(seqs)] <- paste(seq_groups[length(seqs)], "(m)")
  
  # rel_abun_p <- ( (rel_abun_sum - min(rel_abun_sum)) / (max(rel_abun_sum) - min(rel_abun_sum)) ) * (100 - sqrt(max(rel_abun_sum)) - 1) + 1
  rel_abun_p <- rel_abun_sum
  
  data_matrix$seq_plot_pos <- data_matrix$seq
  levels(data_matrix$seq_plot_pos) <- rel_abun_p
  data_matrix$rel_abun <- data_matrix$seq
  levels(data_matrix$rel_abun) <- rel_abun
  print(c("rel_abun", rel_abun))
  print(c("rel_abun_sum", rel_abun_sum))
  print(c("rel_abun_p", rel_abun_p))
  

  gg <- ggplot(data_matrix, aes(x=position, y=as.numeric(as.character(seq_plot_pos)))) +
    # plot horizontal lines with relative abundances
    geom_hline(yintercept = rel_abun_p, 
              size = sqrt(as.numeric(as.character(rel_abun))), 
               color = "grey") +
    # plot vertical lines for mismatches
    geom_point(shape = "|", 
               aes(colour = value, 
                   size=sqrt(as.numeric(as.character(rel_abun))))) +
    # prevent dropping of master sequence despite no points
    scale_y_continuous(limits = c(min(rel_abun_sum)-1, max(rel_abun_sum) + 1), breaks=rel_abun_p, labels = seq_groups) +
    scale_x_continuous(limits = c(1, NA)) +
    labs(x = "Alignment Position", y = element_blank(), title = "Mismatches compared to master", subtitle = data$path) +
    theme_linedraw() +
    theme(axis.text.x = element_text(size = rel(2)), axis.title.x = element_text(size = rel(3))) +
    scale_size_identity()
  
  gg
  # if (inherits(data, "nucleotide")) {
  #   gg + scale_color_manual(name = "Legend",
  #                           breaks = c("A", "C", "G", "T", "-"),
  #                           labels = c("A", "C", "G", "T", "Gap"),
  #                           values = c("A" = "#00BF7D", "C" = "#00B0F6", "G" = "#A3A500", "T" = "#F8766D", "-" = "dark grey"))
  # }
  # else if (inherits(data, "amino acid")) {
  #   # TODO: custom colouring
  #   gg
  # }
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

calc_seq_diff <- function(data, master, method = "svn") {
  # @arg data highlineR Data object to calculate sequence differences
  # @arg master sequence to which others will be compared
  
  data$seq_diff <- matrix(ncol = nchar(data$raw_seq[[1]]$sequence), nrow = length(ls(data$compressed))-1)
  
  master_seq <- strsplit(master, "")[[1]]
  
  row_num <- 1
  row_names <- NULL
  
  if (method == "tvt") {
    for (comp in ls(data$compressed)) { # for each sequence in environment
      if (comp != master) { # ignore master
        comp_seq <- strsplit(comp, "")[[1]]
        
        for (i in 1:length(master_seq)) { # for each positon
          if(master_seq[i] != comp_seq[i]) { # check if different
            if ((identical(sort(c(master_seq[i], comp_seq[i])), sort(c("A", "G")))) || (identical(sort(c(master_seq[i], comp_seq[i])), sort(c("C", "T"))))){
              data$seq_diff[row_num,i] <- "transition"
            }
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
  else if(method == "svn") {
    for (comp in ls(data$compressed)) { # for each sequence in environment
      if (comp != master) { # ignore master
        comp_seq <- strsplit(comp, "")[[1]]
        n <- 1
        
        while ((n + 2) < length(master_seq)) {
          master_codon <- paste(master_seq[n:(n+2)], collapse = "")
          comp_codon <- paste(comp_seq[n:(n+2)], collapse = "")
          if (master_codon != comp_codon) {
            if (identical(aaLookup(master_codon), aaLookup(comp_codon))) {
              data$seq_diff[row_num, (n - 1 + which(strsplit(master_codon, "")[[1]] != strsplit(comp_codon, "")[[1]]))] <- "synonymous"
            }
            else{
              data$seq_diff[row_num,(n - 1 + which(strsplit(master_codon, "")[[1]] != strsplit(comp_codon, "")[[1]]))] <- "non-synonymous"
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
  else{
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