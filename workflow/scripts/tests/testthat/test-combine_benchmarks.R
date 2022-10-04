## This test file is designed to be run from workflow/scripts
## with the command `testthat::test_file("tests/testthat/test-combine_benchmarks.R", reporter = default_reporter())`

source("../../combine_benchmarks.R")

col.types <- cols(
  s = col_double(),
  `h:m:s` = col_time(format = ""),
  max_rss = col_double(),
  max_vms = col_double(),
  max_uss = col_double(),
  max_pss = col_double(),
  io_in = col_double(),
  io_out = col_double(),
  mean_load = col_double(),
  cpu_time = col_double()
)

col.names <- c("s","h:m:s","max_rss","max_vms","max_uss","max_pss","io_in","io_out","mean_load",	"cpu_time", "process", "rule")

test.file1 <- "testthat_resources/benchmarks/21small.tsv"
test.file2 <- "testthat_resources/benchmarks/21small_1_of_3.tsv"

test_that("read.benchmarks properly reads in tsv file", {
    #create expected out df
    vals <- c(5.6752, 5, 113.45, 2586, 107,	108, 80.51, 0.00,	72.63,	0.2, "21small", "benchmarks")
    tempdf <- data.frame(matrix(nrow=1, data=vals))
    tbl <- as_tibble(tempdf)
    tbl$X2 = hms::hms(seconds = 5, minutes = 0)
    exp.df <- tbl %>% mutate_at(vars(-c("X2", "X11", "X12")), as.double)
    colnames(exp.df) <- col.names

    test.out <- read.benchmarks(test.file1)

    expect_identical(test.out, exp.df)
})
