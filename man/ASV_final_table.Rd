% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ASV_final_table}
\alias{ASV_final_table}
\title{ASV table from olive trees endosphere samples.}
\format{
A data frame with 1,065 rows and 32  variables: \describe{
\item{ASV_seqs}{sequences infered with dada2}
\item{Kingdom}{classification} \item{Phylum}{classification}
\item{Class}{classification}
\item{Order}{classification}\item{Family}{classification}
\item{Genus}{classification} \item{ASV_names}{unique names for ASVs with
numbers}\item{L1-1}{From L1-1 to L1-8, samples of olive trees from location
\enumerate{
\item The same for L2 and L3}... }
}
}
\source{
data from our laboratory.
}
\usage{
ASV_final_table
}
\description{
Data frame containing ASV sequences, names, classification and abundance for
olive trees endosphere samples. This ASV table has been cleaned from
eukaryotic sequences and a cut-off has been applied from MOC community data.
}
\keyword{datasets}
