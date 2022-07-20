## Load packages
library(tibble)
library(tidyr)
library(readr)

## Set simulation parameters
debug <- FALSE

if(debug){
  nrep <- 1
}
if(!debug){
  nrep <- 50
}

## Define simulation parameters
pars_mat <- crossing(
  N = c(10, 25, 50),
  tau = c(0.25, 0.35, 0.45, 0.55, 0.65),
  nsites = 4**2,
  nsample_cap = 5,
  nsample_pa = 7,
  ntraps = 3
) %>%
  rowid_to_column("sim") %>%
  crossing(rep = 1:nrep) %>%
  add_column(seed = as.integer(1e6 * runif(nrow(.)))) %>%
  rowid_to_column("task")

write_csv(pars_mat, "./parameter_matrix.csv")
