#' Determine domain isotype
#' @param pfamRes
#' A data frame with pfam augmneted with indel and truncation info (as produced by \code{augment_pfam}).
#' @param fracCutoff
#' The fraction of a protein domain that must be affected before classifying it a truncation or indel.
#' @return
#' The data.frame with the Pfam results now augmented with info about domain domain isotype
#' @examples
#' ### Load pfam data
#' pfamResultFile <- system.file("extdata/pfam_results.txt", package = "pfamAnalyzeR")
#' pfamRes <- read_pfam(pfamResultFile)
#'
#' ### Augment the pfam data
#' pfamRes <- augment_pfam(pfamRes)
#'
#' ### Predict domain isotype
#' pfamRes <- analyse_pfam_isotypes(pfamRes)

analyse_pfam_isotypes <- function(
    pfamRes,
    fracCutoff = 0.1
) {
    pfamRes %>%
        mutate(
            domain_isotype = case_when(
                abs(indel_frac) <=  fracCutoff & truncation_frac <= fracCutoff ~ 'Reference',
                abs(indel_frac)  >  fracCutoff & truncation_frac  > fracCutoff ~ 'Complex',
                indel_frac       >  fracCutoff & truncation_frac <= fracCutoff ~ 'Insertion',
                indel_frac       < -fracCutoff & truncation_frac <= fracCutoff ~ 'Deletion',
                abs(indel_frac) <=  fracCutoff & truncation_frac  > fracCutoff ~ 'Truncation',
                TRUE ~ 'Something went wrong - contact author'
            ),
            domain_isotype_simple =case_when(
                domain_isotype == 'Reference' ~ 'Reference',
                domain_isotype %in% c('Complex','Insertion','Deletion','Truncation') ~ 'Non-reference',
                TRUE ~ 'Something went wrong - contact author'
            )
        )
}
