# Load packages
library(coda)
library(rjags)

print("finished loading libraries")

# source files
source("code/data_generation_funcs.R") # source other functions you have written

# generate data

# number of sites
nsites <- 4**2

# number of individuals in capture history
n <- 10

# number of sampling occasions in capture-recapture
nsample_cap <- 5

# number of sampling occasions in presence-absence
nsample_pa <- 7

# amount of movement for the species
tau <- 0.3

# number of camera traps
ntraps <- 3

# homeranges for the n individuals
homeranges <- create_homeranges(n)

capture_hist <- generate_capture_hist(ntraps, nsample_cap, tau, homeranges)
capture_hist

presence_absence <- generate_presence_absence(nsites, nsample_pa, tau, homeranges)
presence_absence

output_data(presence_absence, capture_hist)