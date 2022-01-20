library(R.utils, quietly=TRUE)
library(stringr, quietly=TRUE)

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

#' Add FastQC fail annotation by sample/lane to existing subject summary data frame
#'
#' @param df data.frame; data frame with subject ID in first column
#' @param fastqc.filename character vector; name of multiqc backend summary
#' table of FastQC output
#' @return data.frame; input df with additional summary information added as new columns
add.fastqc.data <- function(df, fastqc.filename) {
    stopifnot(is.data.frame(df))
    stopifnot(is.character(fastqc.filename))
    stopifnot(file.exists(fastqc.filename))
    fastqc.data <- read.table(fastqc.filename, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE, comment.char = "")
    target.cols <- c("per_base_sequence_quality",
                     "per_base_n_content",
                     "overrepresented_sequences")
    names(target.cols) <- c("Per Base Sequence Quality Failures",
                            "Per Base N Content Failures",
                            "Overrepresented Sequences Failures")
    res <- df
    for (i in seq_len(length(target.cols))) {
        stopifnot(!(names(target.cols)[i] %in% colnames(res)))
        stopifnot(target.cols[i] %in% colnames(fastqc.data))
        res[, names(target.cols)[i]] <- unname(sapply(res[, 1], function(id) {
            target.rows <- stringr::str_detect(fastqc.data[, 1],
                                               paste("^", id, "_S[0-9]+_L[0-9]+_r[12]$", sep = "")) &
                fastqc.data[, target.cols[i]] == "fail"
            if (length(which(target.rows)) == 0) {
                ""
            } else {
                paste(sort(stringr::str_replace(fastqc.data[target.rows, 1],
                                                paste("^", id, "_", sep = ""),
                                                "")), collapse = ", ")
            }
        }))
    }
    res
}

#' Create a table reporting pairs of subjects related
#' to some degree according to Somalier
#'
#' @param output.subjects.filename character vector, name
#' of file containing subject list from `bcftools query -l`
#' @param somalier.relatedness.filename character vector,
#' Somalier paired relatedness report filename
#' @param min.cutoff numeric, lower exclusive bound of
#' Somalier relatedness metric for related pairs
#' @param max.cutoff numeric, upper inclusive bound of
#' Somalier relatedness metric for related pairs
#' @return data.frame, first three columns of Somalier
#' paired relatedness report for subject pairs matching
#' requested relatedness criteria
report.related.subject.pairs <- function(output.subjects.filename, somalier.relatedness.filename, min.cutoff, max.cutoff) {
    stopifnot(is.character(output.subjects.filename))
    stopifnot(file.exists(output.subjects.filename))
    stopifnot(is.character(somalier.relatedness.filename))
    stopifnot(file.exists(somalier.relatedness.filename))
    stopifnot(is.numeric(min.cutoff))
    stopifnot(is.numeric(max.cutoff))
    stopifnot(min.cutoff <= max.cutoff)

    output.subjects <- read.table(output.subjects.filename, header = FALSE, stringsAsFactors = FALSE)[, 1]
    somalier.relatedness <- read.table(somalier.relatedness.filename, header = FALSE, stringsAsFactors = FALSE)[, 1:3]

    somalier.relatedness <- somalier.relatedness[somalier.relatedness[, 1] %in% output.subjects &
                                                 somalier.relatedness[, 2] %in% output.subjects &
                                                 somalier.relatedness[, 3] > min.cutoff &
                                                 somalier.relatedness[, 3] <= max.cutoff, ]
    rownames(somalier.relatedness) <- NULL
    colnames(somalier.relatedness) <- c("Subject 1",
                                        "Subject 2",
                                        "Somalier Relatedness")
    somalier.relatedness
}


#' Pull coverage information from bcftools stats output run on final
#' merged vcf, and use it to populated a column in the summary table
#'
#' @param df data.frame, contains the samples and summary metrics as populated
#' by other functions in this script
#' @param bcftools.stats.filename character vector, filename of text output
#' of bcftools stats
#' @return data.frame, previous summary stats plus coverage as reported by
#' bcftools
add.coverage <- function(df, bcftools.stats.filename) {
    stats.lines <- readLines(bcftools.stats.filename)
    psc.lines <- stats.lines[str_detect(stats.lines, "^PSC")]
    cvg.df <- data.frame(t(data.frame(lapply(str_split(psc.lines, "\t"), function(x) {x[c(3, 10)]}))))
    rownames(cvg.df) <- cvg.df[, 1]
    df[, "Coverage"] <- cvg.df[df[, 1], 2]
    df
}
