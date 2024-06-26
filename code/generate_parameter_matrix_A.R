## Load packages
library(tibble)
library(tidyr)
library(readr)

## Set simulation parameters
debug <- FALSE

## Define simulation parameters
pars_mat <- crossing(
  N = c(10, 50, 100, 200, 500),
  tau = c(0.15, 0.25, 0.35, 0.45, 0.55),
  nsites = c(8, 9, 10)**2,
  nsample_cap = 7,
  nsample_pa = 5,
  ntraps = c(5,15,30),
  p = .75,
  sigma = 1,
  delta = c(1, 4, 8, 16)
) %>%
  rowid_to_column("sim") %>%
  rowid_to_column("task")

write_csv(pars_mat, "./parameter_matrix_A.csv")


