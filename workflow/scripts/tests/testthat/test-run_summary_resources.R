## This test file is designed to be run from workflow/scripts
## with the command `testthat::test_file("tests/testthat/test-run_summary_resources.R", reporter = default_reporter())`

source("../../run_summary_resources.R")

input.subjects <- c("A", "B", "C", "D", "E", "F", "G")
output.subjects.filename <- "testthat_resources/output_subjects.tsv"
exclude.reasons.filename <- "testthat_resources/exclude_reasons.tsv"
somalier.relatedness.filename <- "testthat_resources/somalier.pairs.tsv"
somalier.sex.filename <- "testthat_resources/somalier.samples.tsv"
fastqc.filename <- "testthat_resources/multiqc_fastqc_1.txt"
bcftools.stats.filename <- "testthat_resources/joint_called_stats.out"

test_that("count.rows.in.file works correctly", {
    expect_equal(count.rows.in.file(output.subjects.filename),
                 3)
})

test_that("prepare.subject.tracking.table correctly aligns subjects", {
    expected <- data.frame(Subject = input.subjects,
                           "QC Outcome for This Run" = c("Pass", "low_depth", "Pass", "contamination,low_depth", "No", "Pass", "high_het_hom"),
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

test_that("add.fastqc.data returns lane and read for failing fastq files for specific metrics", {
    input.df <- data.frame(Subject = input.subjects,
			    "QC Outcome for This Run" = c("Pass", "low_depth", "Pass", "contamination,low_depth", "No", "Pass", "high_het_hom"),
			    check.names = FALSE)
    expected <- cbind(input.df,
			"Per Base Sequence Quality Failures" = rep("", nrow(input.df)),
			"Per Base N Content Failures" = rep("", nrow(input.df)),
			"Overrepresented Sequences Failures" = rep("", nrow(input.df)),
			"Rerun Recommendation" = c("Rerun subset of fastqs", "Fail", "Pass", "Fail", "Fail", "Pass", "Fail"))
    expected[1, 3] = "S1_L002_r1, S1_L002_r2"
    expected[1, 4] = "S1_L002_r1"
    expected[2, 5] = "S1_L001_r1"
	observed <- add.fastqc.data(input.df, fastqc.filename)
    expect_identical(observed, expected)
})

test_that("add.coverage returns dataframe with additional column enumerating coverage values per subject", {
    input.df <- data.frame(Subject = input.subjects,
			"QC Outcome for This Run" = c("Pass", "low_depth", "Pass", "contamination,low_depth", "No", "Pass", "high_het_hom"),
			check.names = FALSE)
    expected <- input.df
	expected$Coverage <- c(18.5, 22.4, 23.6, 29.1, 28.7, 25.2, 19.1)
    observed <- add.coverage(input.df, bcftools.stats.filename)
	expect_identical(observed, expected)
})

test_that("report.sex.discordances returns dataframe with subject, inferred, and self-reported sex, if discordant", {
	observed <- report.sex.discordances(output.subjects.filename, somalier.sex.filename)
	expected <- data.frame("Subject" = c("C", "F"), "Inferred Sex" = c(1, 2), "Self-reported Sex" = c("-9", "male"), check.names = FALSE)
	expect_equal(observed, expected)
})
