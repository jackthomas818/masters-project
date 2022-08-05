# Load packages
library(coda)
library(rjags)
library(tictoc)

# R script for running a simulation study. Uses a parameter matrix
# created by generate_parameter_matrix.R. Takes the task id as a command
# line argument,which corresponds to a row of the parameter matrix.
# JAGS code taken from Blanc et. al (2014). This simulation study
# intends to test the model suggested by Blanc et. al (2014). Outputs
# model_statistics_* text files with the summary statistics for the coda
# samples.


tic()

print("finished loading libraries")

# source files
source("./data_generation_funcs.R")

args <- commandArgs(trailingOnly = TRUE)

task_id_str <- args[1]
task_id <- strtoi(args[1])

# number of repetitions for each simulation
reps <- args[2]

print(paste0("task id: ", task_id_str))

# parameter matrix from generate_parameter_matrix.R
params_matrix <- read.csv("./parameter_matrix.csv")

print("finished loading parameter matrix")

params <- params_matrix[task_id, ]

# generate data

# number of sites
nsites <- as.integer(params["nsites"])

# number of individuals in capture history
n_total <- as.integer(params["N"])

# number of sampling occasions in capture-recapture
nsample_cap <- as.integer(params["nsample_cap"])

# number of sampling occasions in presence-absence
nsample_pa <- as.integer(params["nsample_pa"])

# amount of movement for the species
tau <- as.double(params["tau"])

# number of camera traps
ntraps <- as.integer(params["ntraps"])

for (rep in 1:reps) {
  # setting seed
  seed <- as.integer(1e7 * runif(1))
  set.seed(seed)

  print(paste0("seed: ", seed))

  # homeranges for the n individuals
  homeranges <- create_homeranges(n_total)

  capture_hist <- generate_capture_hist(ntraps, nsample_cap, tau, homeranges)
  capture_hist

  presence_absence <- generate_presence_absence(nsites, nsample_pa, tau, homeranges)
  presence_absence

  output_data(presence_absence, capture_hist, paste0(task_id_str,"-",rep))

  # Supplementary Code by Blanc et al. (2014)

  #---------------------- patch occupancy data ----------------------------#
  # Presence-absence data structure: rows = sites; columns = occasions
  y <- read.csv(paste0("./presence-absence-data-", task_id_str, ".csv"), header = FALSE) # load presence-absence data
  nsites <- dim(y)[1]
  nsurvs <- dim(y)[2]
  #--------------------- capture-recapture data -----------------------------#
  # Capture-recapture data structure : rows = individuals; columns = capture occasions
  mydata <- read.table(paste0("./capture-recapture-data-", task_id_str, ".txt"), header = FALSE) # load capture-recapture data
  extra <- 250 # define large number of extra individual capture histories
  n <- nrow(mydata) # number of observed individuals
  M <- extra + n
  xn <- rowSums(mydata)
  x <- c(xn, rep(0, extra))
  k <- ncol(mydata)
  zerouse <- 0

  #-------------------- Jags stuff -----------------------#
  # List of data
  mydatax <- list(x = x, M = M, n = n, k = k, zerouse = zerouse, y = y, nsites = nsites, nsurvs = nsurvs)

  # Parameters monitored
  parameters <- c("N", "ppo", "mumup0", "psi0", "mupsi", "mean.p0", "sdeps", "taueps", "betapsi", "etap")


  # Initial values
  init1 <- list(mupsi = runif(1), truocc = as.numeric(apply(y, 1, sum) > 0), xi = rnorm(1, sd = 2), N = sample(n:M, 1), mumup = rnorm(1), taup = rlnorm(1))
  init2 <- list(mupsi = runif(1), truocc = as.numeric(apply(y, 1, sum) > 0), xi = rnorm(1, sd = 2), N = sample(n:M, 1), mumup = rnorm(1), taup = rlnorm(1))
  init3 <- list(mupsi = runif(1), truocc = as.numeric(apply(y, 1, sum) > 0), xi = rnorm(1, sd = 2), N = sample(n:M, 1), mumup = rnorm(1), taup = rlnorm(1))
  inits <- list(init1, init2, init3)

  #-------------------- Call jags from R -----------------------#
  jmodel <- jags.model("CRandPO_HET.txt", mydatax, inits, n.chains = 3, n.adapt = 2500)
  jsample <- coda.samples(jmodel, parameters, n.iter = 15000, thin = 1)
  save(jsample, file = paste0("coda_samples_", task_id))

  # Get summary statistics
  s <- summary(jsample)

  all_table <- cbind(s$statistics, s$quantiles)

  write.table(all_table, file = paste0("model_statistics_", task_id_str, "_", toString(rep), "_", toString(seed), ".txt"))
}


print(tictoc::toc())
