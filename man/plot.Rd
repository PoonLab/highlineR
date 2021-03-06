% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot}

\alias{plot}
\alias{plot.session}
\alias{plot.Data}

\title{highlightR plotting}
\usage{
\method{plot}{session}(x, mode = "mismatch", master, sort_by, rf = 1,
  use_sample = T, ...)

\method{plot}{Data}(x, session_plot = F, mode = "mismatch",
  master = x$master, sort_by = "similarity", rf = 1,
  use_sample = T, ...)

}
\arguments{
\item{x}{Data object or session object containing Data objects to plot.}

\item{mode}{A character string representing the desired mutation annotation for plotting. Options: "mismatch" (default), "svn" (Synonymous versus Non-Synonymous), "tvt" (Transition versus Transversion).}

\item{master}{A character string representing the sequence to which other sequences should be compared. By default, the most abundant sequence is selected.}

\item{sort_by}{A character string representing how the sequences should be ordered in a plot. Options: "similarity", "frequency".}

\item{rf}{An integer specifying which reading frame should be used when determining Synonymous vs Non-Synonymous mutations. Options: 1 (default), 2, 3.}

\item{use_sample}{A logical value specifying which environment should be plotted. If \code{True}, then the \code{sample} environment of the Data objects is plotted. If \code{False}, the complete \code{compressed} environment is plotted.}

\item{seq_diff}{matrix of compositional differences between sequences.}

\item{seq_order}{string vector of sequences present in compressed environment ordered as desired for plotting.}
}
\value{
Plots NGS data
}
\description{
Plots highlineR Data objects.
}
\keyword{internal}
