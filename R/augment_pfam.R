#' Augment pfam domains with truncation/indel calculations
#' @param pfamRes
#' A data frame with pfam results as produced by \code{read_pfam}.
#' @param useAllignmentCoordinates
#' A logic indicating whether to use the coordinates of the aligned region (if set to TRUE) or the envelope region (if set to FALSE). Defaults to TRUE.
#' @return
#' The data.frame with the Pfam results now augmented with info on trunkation and indel sizes
#' @examples
#' ### Load pfam data
#' pfamResultFile <- system.file("extdata/pfam_results.txt", package = "pfamAnalyzeR")
#' pfamRes <- read_pfam(pfamResultFile)
#'
#' ### Augment the pfam data
#' pfamRes <- augment_pfam(pfamRes)
#'

augment_pfam <- function(
    pfamRes
) {
    ### Calculate truncation
    if(TRUE) {
        startMissing <- pfamRes$hmm_start - 1
        endMissing   <- pfamRes$hmm_length - pfamRes$hmm_end

        startFracMissing <- round( startMissing / pfamRes$hmm_length, digits = 3)
        endFracMissing   <- round( endMissing   / pfamRes$hmm_length, digits = 3)

        pfamRes$truncation_size <- base::pmax(endMissing      , startMissing )
        pfamRes$truncation_frac <- base::pmax(startFracMissing, endFracMissing)
    }

    ### Calculate indel
    if(TRUE) {
        pfamRes <-
            pfamRes  %>%
            tibble::as_tibble() %>%
            dplyr::mutate(
                hmm_alligned_length = hmm_end - hmm_start + 1,
                seq_alligned_length = alignment_end - alignment_start + 1,
                indel_size = seq_alligned_length - hmm_alligned_length,
                indel_frac = round( indel_size / hmm_alligned_length, digits = 3)
            ) %>%
            dplyr::select(-hmm_alligned_length,-seq_alligned_length)

    }
    # indel_size Postive : Sequence analyzed is longer  aka sequence insertion
    # indel_size Negative: Sequence analyzed is shorter aka sequence deletion

    return(pfamRes)
}
