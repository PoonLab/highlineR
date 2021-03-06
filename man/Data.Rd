% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import.R
\name{Data}
\alias{Data}
\title{Constructor function for S3 class representing an NGS Data object.}
\usage{
Data(path, datatype = utils::tail(strsplit(path, "\\\\.")[[1]], n = 1),
  seqtype = "nucleotide")
}
\arguments{
\item{path}{A character string representing the absolute path to sequence file.}

\item{datatype}{A character string representing the file type. Default value is file extension (last of all period (".") separated extensions).}

\item{seqtype}{A character string representing the sequence type. Options: "nucleotide" (default) or "amino acid".}
}
\value{
Returns a constructed S3 object which is a child object of a parent session object, inherits from the following classes:
\itemize{
\item \code{datatype}: The file type of NGS file represeted by the object.
\item \code{seqtype}: The sequence type of the sequences in the NGS file represeted by the object.
\item \code{Data}.
\item \code{environment}.
}
and contains the following elements:
\itemize{
\item \code{path}: The absolute path to the NGS file represented by the object.
\item \code{raw_seq}: A list that will contain the parsed sequences of the NGS file.
\item \code{compressed}: An S3 object that inherits from the base \code{environment} class and will contain the compressed sequences of the NGS file.
\item \code{sample}: An S3 object that inherits from the base \code{environment} class and will contain randomly sampled sequences from the the compressed structure.
\item \code{master}: A character string representing the master sequence to use for plotting.
\item \code{seq_diff}: A character matrix of sequence differences in variant sequences compared to the master sequence.
}
}
\description{
\code{Data} returns constructed Data object.
}
