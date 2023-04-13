#' Read Pfam file into R
#' @description
#' Read Pfam result file file into R. Supports both result files from local and web-server
#' @param path
#' A string indicating the full path to the Pfam result file
#' @details
#' The pfam webserver can be found at \url{https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan}.
#' @return
#' A data.frame with the Pfam results
#' @examples
#' pfamResultFile <- system.file("extdata/pfam_results.txt", package = "pfamAnalyzeR")
#' pfamRes <- read_pfam(pfamResultFile)

read_pfam <- function(
    path
) {

    ### Figure out file type
    if(TRUE) {
        fileVector <- base::readLines(con = path, n = 3)

        if (! length(fileVector) ) {
            stop('The file pointed to by \'path\' is empty')
        }

        tapsAnnotated <- any(stringr::str_detect(fileVector, '\\t'))
    }
    # tapsAnnotated

    ### Figure out header and lines to skip
    if(TRUE) {
        ### Test whether headers are included
        temp <-
            utils::read.table(
                file = path,
                stringsAsFactors = FALSE,
                fill = TRUE,
                header = FALSE,
                nrows = 1
            )
        headerIncluded <- grepl('^<seq|^seq', temp[1, 1])

        ### Figure out number of lines to skip
        if(   headerIncluded) {
            tmp2 <-
                base::readLines(
                    con = path,
                    n = 50
                )
            skipLine <- which(grepl('^<seq|^seq', tmp2))[1]

            temp3 <-
                utils::read.table(
                    file = path,
                    stringsAsFactors = FALSE,
                    fill = TRUE,
                    header = FALSE,
                    nrows = 1,
                    skip = skipLine,
                )
        }
        if( ! headerIncluded) {
            skipLine <- 0
        }
    }
    # skipLine

    ### Read in file using the file type info extracted above
    if(TRUE) {
        if( tapsAnnotated ) {
            localPfamRes <- utils::read.table(
                file = path,
                stringsAsFactors = FALSE,
                fill = TRUE,
                header = FALSE,
                col.names = seq_len(16),
                sep = '\t',
                skip = skipLine
            )
        } else {
            localPfamRes <- utils::read.table(
                file = path,
                stringsAsFactors = FALSE,
                fill = TRUE,
                header = FALSE,
                col.names = seq_len(16),
                sep = '', # one or more white spaces
                skip = skipLine
            )
        }

    }

    ### Test for problems
    if(TRUE) {
        ### Determine if active sites are predicted
        activeResidueIncluded <- any(
            stringr::str_detect(
                string = localPfamRes[,ncol(localPfamRes)],
                pattern = 'predicted_active_site'
            ),
            na.rm = TRUE
        )

        ### Test shifted values
        testShiftedValues <- function(aDF) {
            try1 <- unique(c(
                which(is.na( aDF[,14]     )),              # via "significant" column (will either be NA or a clan indication)
                which(       aDF[,14] != 1 ),              # via "significant" column (will either be NA or a clan indication)
                which( ! stringr::str_detect(aDF[,6], '^PF|^PB') )  # via pfam_hmm id which should start with PF or PB
            ))

            return(try1)
        }
        try1problems <- testShiftedValues(localPfamRes)
        any(try1problems)
    }
    # FALSE

    ### Try with fwf
    if( any(try1problems) ) {
        ### Fix the mistake in the fixed width file
        suppressWarnings(
            fileVector <- readLines(path, encoding = 'UTF-8')
        )
        fileVector <- gsub('Coiled-coil',' Coiled', fileVector)

        ### FOr read_fwf to recognice it as raw text it must have new lines included
        fileVector <- stringr::str_c(fileVector, '\n')

        suppressMessages(
            suppressWarnings(
                localPfamRes2 <- readr::read_fwf(
                    file = fileVector,
                    col_positions = readr::fwf_empty(
                        file = fileVector,
                        col_names = paste0(
                            'X',
                            seq_len(15+ activeResidueIncluded + 1)
                        ),
                        skip = skipLine,
                        comment = '#',
                        n = 1000L
                    ),
                    skip = skipLine,
                    comment = '#'
                )
            )
        )
        localPfamRes2 <- as.data.frame(localPfamRes2)

        ### Test shifted values via "significant" column (will either be NA or a clan indication)
        try2problems <- testShiftedValues(localPfamRes2)

        ### overwrite if fewer problems
        if( length(try2problems) < length(try1problems) ) {
            localPfamRes <- localPfamRes2
        }
    }

    ### Set col names
    if(TRUE) {
        ### Define Colnames for pfam file types
        if(TRUE) {
            oldColnames <-
                c(
                    'seq_id',
                    'alignment_start',
                    'alignment_end',
                    'envelope_start',
                    'envelope_end',
                    'hmm_acc',
                    'hmm_name',
                    'type', # 8
                    'hmm_start',
                    'hmm_end',
                    'hmm_length',
                    'bit_score',
                    'e_value',
                    'significant', # 14
                    'clan',
                    'residue'
                )
            newColNames <-
                c(
                    'seq_id',
                    'alignment_start',
                    'alignment_end',
                    'envelope_start',
                    'envelope_end',
                    'hmm_acc',
                    'hmm_name',
                    'hmm_start', # 8
                    'hmm_end',
                    'hmm_length',
                    'bit_score',
                    'Individual_E_value',
                    'Conditional_E_value',
                    'significant', # 14
                    'outcompeted',
                    'clan'
                )
        }

        ### Define style via clan column
        oldStyle <- any(grepl('^CL', localPfamRes[,15]))
        if(ncol(localPfamRes) == 16 ) {
            newStyle <- any(grepl('^CL', localPfamRes[,16]))
        } else {
            newStyle <- FALSE
        }

        ### Old style
        if ( oldStyle & !newStyle ) {
            colnames(localPfamRes) <- oldColnames[seq_len(ncol(localPfamRes))]

            if( 'residue' %in% colnames(localPfamRes) ) {
                localPfamRes$residue[which(localPfamRes$residue == '')] <-
                    NA
            }
        }
    }

    return(localPfamRes)
}
