library(R.utils, quietly=TRUE)

#' Determine the number of rows in a file
#'
#' @param filename character vector; name of file to probe
#' @return integer count of rows
#'
#' Note that R.utils masks an enormous number of functions,
#' including some used by testthat, so this may eventually
#' cause problems downstream.
count.rows.in.file <- function(filename) {
    stopifnot(is.vector(filename, mode = "character"),
              length(filename) == 1)
    stopifnot(file.exists(filename))
    R.utils::countLines(filename)[[1]]
}


#' Prepare a table reporting subject resolution
#' through the entire pipeline.
#'
#' @param input.subjects character vector; input
#' subject IDs from manifest, without duplicates
#' @param output.subjects.filename character vector;
#' file listing subjects making it through the pipeline,
#' one subject ID per file
#' @return data.frame, prepared table ready for knitr
#'
prepare.subject.tracking.table <- function(input.subjects, output.subjects.filename, exclude.reasons.filename) {
    stopifnot(is.vector(input.subjects, mode = "character"))
    stopifnot(is.vector(output.subjects.filename, mode = "character"),
              length(output.subjects.filename) == 1)
    stopifnot(file.exists(output.subjects.filename))
    stopifnot(file.exists(exclude.reasons.filename))
    exclude.reasons <- read.table(exclude.reasons.filename, header = FALSE, sep = "\t", row.names = 1)
    output.subjects <- read.table(output.subjects.filename, header = FALSE, sep = "\t")[, 1]
    result <- rep("Pass", length(input.subjects))
    result[!(input.subjects %in% output.subjects)] <- "No"
    result[input.subjects %in% rownames(exclude.reasons)] <- exclude.reasons[input.subjects[input.subjects %in% rownames(exclude.reasons)], 1]
    df <- data.frame("Subject" = input.subjects,
                     "Final QC Outcome" = result,
                     check.names = FALSE)
    df
}
