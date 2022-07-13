## Load packages
library(tidyverse)

## Set simulation parameters
debug <- TRUE

if(debug){
  nrep <- 2
}
if(!debug){
  nrep <- 50
}

## Define simulation parameters
pars_mat <- crossing(
  N = c(10, 25, 50, 100),
  tau = c(0.1, 0.2, 0.3, 0.4, 0.5),
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
