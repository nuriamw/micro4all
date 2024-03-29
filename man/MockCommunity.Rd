% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mockcommunity.R
\name{MockCommunity}
\alias{MockCommunity}
\title{Apply a cut-off based on Mock Community sequencing results}
\usage{
MockCommunity(data, MOCK_composition, ASV_column, choose.first = FALSE)
}
\arguments{
\item{data}{a data frame, \strong{rows corresponding to OTUs} and
\strong{columns to classification}, followed by sample composition. MOCK
community samples should include the word \emph{MOCK} in their names.}

\item{MOCK_composition}{a data frame containing \strong{classification for the
MOCK community}, usually provided by the company. The first column is the
one to be used by the function, corresponding to genus level.}

\item{ASV_column}{name of the column in \code{data} where ASV or OTU names are
stored.}

\item{choose.first}{whether to choose first spurious ASV found to calculate
the percentage or not. Default is FALSE.}
}
\value{
Returns a data frame with the same structure as \code{data}, but with
a \strong{cut-off applied based on MOCK community information} and with MOCK
samples removed.
}
\description{
This function reads an OTU table in a specific format and remove
those ASV which represent less than a calculated percentage of sequences.
This percentage is calculated based on MOCK community members (more
information in details).
}
\details{
This function is supposed to be used in the following situation.
Imagine a  metabarcoding experiment for the study of microbial communities.
In the same run, some samples of a commercial MOCK community are analysed.
Then, they are processed at bioinformatic level, alongside the main
experiment samples. The MOCK community results are analyzed through this
function, comparing its composition with the real composition (known design
or provided by the company). When a microorganism is found that should not
be in the MOCK community, it is attributed to sequencing errors. Based on
the abundance of this microorganism (ASV or OTU), the function calculates
the \strong{percentage of sequences} that it represents to the total of the
MOCK community samples. This percentage is used as a cut-off for the main
experiment samples, where ASVs or OTUs representing less than this
percentage are removed.
}
\examples{

ASV_table_cleaned <- MockCommunity(ASV_table_classified_raw, mock_composition, "ASV_names",
choose.first=TRUE)

}
