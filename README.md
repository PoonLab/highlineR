
<!-- README.md is generated from README.Rmd. Please edit that file -->

# highlineR: a tool for visualizing NGS datasets in R

[Highlighter](https://www.hiv.lanl.gov/content/sequence/HIGHLIGHT/highlighter_top.html)
is a web-based tool maintained by the Los Alamos National laboratory
that makes it easy for investigators to quickly scan their sequence
alignments for compositional differences. However, it is inconvenient
for users to use Highlighter to process a large number of alignment
files as the tool only accepts a single file as input and limits users
to 500 sequences per file. This limit on the number of sequences can
also be problematic for users working with alignments derived from
next-generation sequencing platforms. Furthermore, the tool does not
visualize the frequencies of different sequence variants in a dataset.

The *highlineR* package provides users with easy-to-use functions for
batch-processing of sequence alignments, including NGS data, to generate
plots similar to Highlighter. In addition, *highlineR* uses the visual
aid of varying line widths to indicate variant frequencies in the data
(see *Usage*).

## Installation

The simplest method to install *highlineR* is to use the R *devtools*
package:

``` r
# install.packages("devtools")  # if not already installed
require(devtools)
devtools::install_github("PoonLab/highlineR")
```

## Usage

``` r
require(highlineR)
#> Loading required package: highlineR
files <- Sys.glob('~/git/highlineR/inst/extdata/*.fasta')
highline(files[1:2])  # render the first two alignments
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" style="display: block; margin: auto;" />

Note that for the above example, the files were loaded from a developer
directory. To load these same files from your installed package as a
user, you’d have to replace the following line:

``` r
# files <- Sys.glob('~/git/highlineR/inst/extdata/*.fasta')
files <- Sys.glob(paste0(system.file(package='highlineR'), '/extdata/*.fasta'))
```

## Explanation

A *multiple sequence alignment* (MSA) is a hypothesis about how residues
(nucleotides or amino acids) in two or more sequences were derived from
residues in a common ancestral sequence. By convention, each sequence is
arranged horizontally in rows to be read from left to right, and
evolutionarily-related (homologous) residues are arranged vertically in
a column. It is possible for two or more sequences to be exactly
identical. When this occurs, we refer to the shared sequence as the
“sequence variant”, and the number of times this variant appears in
the data as its “count”. (Note there is no established terminology for
these features.) It is common practice to reduce a sequence alignment to
the unique variants, especially when certain variants predominate the
sample population. By doing so, however, we lose information on the
relative abundance of variants.

#### Comparing to a reference sequence

If an MSA comprises a large number of sequences, it may become difficult
to visualize the composition of the alignment. *Highlighter* plots were
developed to reduce the information in the MSA by marking the location
of residues that are different from the reference sequence, which has
been selected by the user or the program. We refer to this reference
sequence the *master sequence*. By default, the `highline` function
selects the most common sequence variant to be the master.
