library(R.utils)

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
prepare.subject.tracking.table <- function(input.subjects, output.subjects.filename) {
    stopifnot(is.vector(input.subjects, mode = "character"))
    stopifnot(is.vector(output.subjects.filename, mode = "character"),
              length(output.subjects.filename) == 1)
    stopifnot(file.exists(output.subjects.filename))
    output.subjects <- read.table(output.subjects.filename, header = FALSE, sep = "\t")[, 1]

    df <- data.frame(Subject = input.subjects,
                     "In Final VCF" = ifelse(input.subjects %in% output.subjects, "yes", "no"),
                     check.names = FALSE)
    df
}
