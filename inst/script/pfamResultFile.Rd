\name{pfamResultFile}
\alias{pfamResultFile}
\docType{data}

\title{
Example of Pfam webserver output
}

\description{
The output from the pfam web-server when running it on 135 protein sequences.
}

\details{
The 135 protein sequences were uploaded to the pfam webserver \url{https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan} as an amino acid fasta file. An email was provided and the results here are copy/pasted from the email they send (as the website at the time of writing have no way to download all results).
}

\source{
The amino acid sequnces analyzed via the pfam server are from the example data in the IsoformSwitchAnalyzeR R package. More information on the data and how it was generated can be found here \code{?IsoformSwitchAnalyzeR::exampleSwitchListIntermediary}.

This data is only included to provide examples for usage of function. As it is modified to illustrate the package it should not be considered real and no biological conclusions should be made from it.
}

\examples{
pfamResultFile <- system.file(
    "extdata/pfam_results.txt",
    package = "pfamAnalyzeR"
)
pfamRes <- read_pfam(pfamResultFile)
}
