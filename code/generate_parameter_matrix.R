## Load packages
library(tibble)
library(tidyr)
library(readr)

## Set simulation parameters
debug <- FALSE

## Define simulation parameters
pars_mat <- crossing(
  N = c(50,200),
  tau = c(0.15,0.20,0.25,0.30,0.35,0.4,0.45,0.5,0.55),
  nsites = c(7**2,8**2,9**2,10**2),
  nsample_cap = c(8,10,12,14),
  nsample_pa = c(8,10,12,14),
  ntraps = c(3,4,5,6)
) %>%
  rowid_to_column("sim") %>%
  rowid_to_column("task")

write_csv(pars_mat, "./parameter_matrix.csv")
