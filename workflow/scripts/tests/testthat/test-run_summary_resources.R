## This test file is designed to be run from workflow/scripts
## with the command `testthat::test_file("tests/testthat/test-run_summary_resources.R", reporter = default_reporter())`

source("../../run_summary_resources.R")

input.subjects <- c("A", "B", "C", "D", "E", "F", "G")
output.subjects.filename <- "testthat_resources/output_subjects.tsv"
exclude.reasons.filename <- "testthat_resources/exclude_reasons.tsv"

test_that("count.rows.in.file works correctly", {
    expect_equal(count.rows.in.file(output.subjects.filename),
                 3)
})

test_that("prepare.subject.tracking.table correctly aligns subjects", {
    expected <- data.frame(Subject = input.subjects,
                           "Final QC Outcome" = c("Pass", "low_depth", "Pass", "contamination,low_depth", "No", "Pass", "high_het_hom"),
                           check.names = FALSE)
    expect_identical(prepare.subject.tracking.table(input.subjects, output.subjects.filename, exclude.reasons.filename),
                     expected)
})
