#' Read in and analyze pfam domains isotypes
#' @param path
#' A string indicating the full path to the Pfam result file
#' @param fracCutoff
#' The fraction of a protein domain that must be affected before classifying it a truncation or indel.
#' @return
#' The data.frame with the Pfam results now augmented with info about domain structural variation
#' @examples
#' ### Predict domain isotypes in pfam results
#' pfamResultFile <- system.file("extdata/pfam_results.txt", package = "pfamAnalyzeR")
#' pfamRes <- pfamAnalyzeR(pfamResultFile)

pfamAnalyzeR <- function(
    path,
    fracCutoff = 0.1
) {
    read_pfam(path) %>%
        augment_pfam() %>%
        analyse_pfam_isotypes(
            pfamRes = .,
            fracCutoff = fracCutoff
        ) %>%
        return()
}
