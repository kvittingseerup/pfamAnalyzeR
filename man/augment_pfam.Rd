% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/augment_pfam.R
\name{augment_pfam}
\alias{augment_pfam}
\title{Augment pfam domains with truncation/indel calculations}
\usage{
augment_pfam(pfamRes)
}
\arguments{
\item{pfamRes}{A data frame with pfam results as produced by \code{read_pfam}.}

\item{useAllignmentCoordinates}{A logic indicating whether to use the coordinates of the aligned region (if set to TRUE) or the envelope region (if set to FALSE). Defaults to TRUE.}
}
\value{
The data.frame with the Pfam results now augmented with info on trunkation and indel sizes
}
\description{
Augment pfam domains with truncation/indel calculations
}
\examples{
### Load pfam data
pfamResultFile <- system.file("extdata/pfam_results.txt", package = "pfamAnalyzeR")
pfamRes <- read_pfam(pfamResultFile)

### Augment the pfam data
pfamRes <- augment_pfam(pfamRes)

}
