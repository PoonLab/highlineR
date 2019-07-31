\name{highline}
\alias{highline}

\title{highlineR plotting}

\description{
  \code{highline} generates \href{https://www.hiv.lanl.gov/content/sequence/HIGHLIGHT/highlighter_top.html}{Highlighter}-style visualizations from one or more sequence alignments.
  These plots use colour marks to annotate the differences of each sequence in the alignment from a reference sequence. 
}

\usage{
highline(file, datatype='fasta', seqtype='nucleotide', unique=FALSE, 
         rf=1, mode='mismatch', sort_by='similarity', quiet=TRUE)
}

\arguments{
  \item{file}{a character vector containing one or more relative or absolute paths to the files containing sequence alignment(s).}
  \item{datatype}{specifies the file format - one of 'fasta', 'fastq' or 'csv'.}
  \item{seqtype}{specifies the type of sequence - 'nucleotide' or 'amino acid'}
  \item{unique}{if FALSE, gather identical sequences into unique variants}
}