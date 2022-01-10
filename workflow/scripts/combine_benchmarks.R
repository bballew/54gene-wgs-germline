# load in required libraries
library(tidyverse)

# read tsv input and append rule name and process
read_benchmarks <- function(filename) {
  read_tsv(filename) %>%
    mutate( 
      process = gsub("*.tsv","",basename(filename)),
      rule =  basename(dirname(filename))
      )
}

# run function for all specified snakemake inputs and output tsv
dataset <- snakemake@input[["tsv"]] %>%
  lapply(read_benchmarks) %>%
  bind_rows()
head(dataset)
write_tsv(dataset, snakemake@output[['benchmarks']])