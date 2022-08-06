## Load packages
library(tibble)
library(tidyr)
library(readr)

## Set simulation parameters
debug <- FALSE

## Define simulation parameters
pars_mat <- crossing(
  N = c(10, 25, 50,75,100,125,150,175,200),
  tau = c(0.15,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55),
  nsites = 4**2,
  nsample_cap = 7,
  nsample_pa = 5,
  ntraps = c(3,4,5)
) %>%
  rowid_to_column("sim") %>%
  rowid_to_column("task")

write_csv(pars_mat, "./parameter_matrix.csv")
