## This test file is designed to be run from workflow/scripts
## with the command `testthat::test_file("tests/testthat/test-run_summary_resources.R", reporter = default_reporter())`

source("../../run_summary_resources.R")

input.subjects <- c("A", "B", "C", "D", "E", "F", "G")
output.subjects.filename <- "testthat_resources/output_subjects.tsv"
exclude.reasons.filename <- "testthat_resources/exclude_reasons.tsv"
somalier.relatedness.filename <- "testthat_resources/somalier.pairs.tsv"

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

test_that("report.related.subject.pairs returns requested subject set", {
    expected <- data.frame("Subject 1" = c("A", "C"),
                           "Subject 2" = c("C", "F"),
                           "Somalier Relatedness" = c(0.4, 0.99),
                           check.names = FALSE)
    observed <- report.related.subject.pairs(output.subjects.filename,
                                             somalier.relatedness.filename,
                                             0.35,
                                             1.0)
    expect_identical(observed, expected)
})

test_that("report.related.subject.pairs returns empty table when range includes no pairs", {
    expected <- data.frame("Subject 1" = c("a"),
                           "Subject 2" = c("b"),
                           "Somalier Relatedness" = c(1.1),
                           check.names = FALSE)
    expected <- expected[-1, ]
    observed <- report.related.subject.pairs(output.subjects.filename,
                                             somalier.relatedness.filename,
                                             0.999,
                                             1.0)
    expect_identical(observed, expected)
})

test_that("report.related.subject.pairs treats min as exclusive and max as inclusive", {
    expected <- data.frame("Subject 1" = c("C"),
                           "Subject 2" = c("F"),
                           "Somalier Relatedness" = c(0.99),
                           check.names = FALSE)
    observed <- report.related.subject.pairs(output.subjects.filename,
                                             somalier.relatedness.filename,
                                             0.4,
                                             0.99)
    expect_identical(observed, expected)
})
