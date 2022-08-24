# Load packages
library(readr)
library(tictoc)
library(data.table)

options(scipen = 999)


#' Calculate Credible Interval Proportion
#'
#' Function that calculates how many credible intervals
#' contain the true value
#'
#' @param lower the array of lower bounds
#' @param upper the array of upper bounds
#' @param actual the true value to be covered
#'
#' @return
#' @export
#'
#' @examples
cred_prop <- function(lower, upper, actual) {
  cbind(lower <= actual & upper >= actual)
}


tic()

print("finished loading libraries")

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# directory of simulation study
directory <- args[1]

# number of repetitions per simulation
reps <- strtoi(args[2])

# directory <- "SimStudy9"
# reps <- 20

# parameter matrix from generate_parameter_matrix.R
params_matrix <- read.csv(paste0(directory, "/parameter_matrix.csv"))

# Getting model statistics
# construct summary statistics table
N_mean <- array(data = NA, dim = nrow(params_matrix) * reps)
N_lower_95 <- array(data = NA, dim = nrow(params_matrix) * reps)
N_upper_95 <- array(data = NA, dim = nrow(params_matrix) * reps)
n_captured <- array(data = NA, dim = nrow(params_matrix) * reps)
N_sd <- array(data = NA, dim = nrow(params_matrix) * reps)
N_naive_se <- array(data = NA, dim = nrow(params_matrix) * reps)
pa_detects_total <- array(data = NA, dim = nrow(params_matrix) * reps)
sim <- array(data = NA, dim = nrow(params_matrix) * reps)
tau <- array(data = NA, dim = nrow(params_matrix) * reps)
nsites <- array(data = NA, dim = nrow(params_matrix) * reps)
ntraps <- array(data = NA, dim = nrow(params_matrix) * reps)
N_actual <- array(data = NA, dim = nrow(params_matrix) * reps)
seed <- array(data = NA, dim = nrow(params_matrix) * reps)

count <- 1
for (task_id in 1:nrow(params_matrix)) {
  print(paste0(task_id, " out of ", nrow(params_matrix)))
  for (rep in 1:reps) {
    model_file <- paste0(directory, "/model_statistics_", task_id, "_", rep, ".txt")
    if (file.exists(model_file)) {
      model_data <- tryCatch(read.table(model_file), error = function(e) NULL)
    } else {
      count <- count + 1
      next
    }

    cr_file <- paste0(directory, "/capture-recapture-data-", task_id, "-", rep, ".txt")
    if (file.exists(cr_file)) {
      cr_data <- tryCatch(read.table(cr_file), error = function(e) NULL)
    } else {
      count <- count + 1
      next
    }

    pa_file <- paste0(directory, "/presence-absence-data-", task_id, "-", rep, ".csv")
    if (file.exists(pa_file)) {
      pa_data <- tryCatch(read.csv(pa_file, header = FALSE), error = function(e) NULL)
    } else {
      count <- count + 1
      next
    }

    # count number of presence detections
    pa_detects_total[count] <- sum(rowSums(pa_data))

    N_mean[count] <- unlist(model_data["Mean"])[1]
    N_lower_95[count] <- unlist(model_data["X2.5."])[1]
    N_upper_95[count] <- unlist(model_data["X97.5."])[1]
    N_sd[count] <- unlist(model_data["SD"])[1]
    N_naive_se[count] <- unlist(model_data["Naive.SE"])[1]
    seed[count] <- unlist(model_data["seed"])[1]

    if (is.null(cr_data)) {
      n_captured[count] <- 0
    } else {
      n_captured[count] <- nrow(cr_data)
    }

    sim[count] <- task_id
    tau[count] <- unlist(params_matrix["tau"])[task_id]
    nsites[count] <- unlist(params_matrix["nsites"])[task_id]
    ntraps[count] <- unlist(params_matrix["ntraps"])[task_id]
    N_actual[count] <- unlist(params_matrix["N"])[task_id]

    count <- count + 1
  }
}

cred_contains_actual_N <- array(cred_prop(N_lower_95, N_upper_95, N_actual))

# calculate bias
N_bias <- array((N_mean - N_actual) / N_actual)

all_data <- data.frame(cbind(sim, N_actual, N_mean, N_bias, n_captured, tau, nsites, ntraps, pa_detects_total, N_lower_95, N_upper_95, cred_contains_actual_N, N_sd, N_naive_se, seed))

# Calculate means by simulation number
setDT(all_data)

summary_stats <- all_data[, list(
  N_actual = mean(N_actual),
  n_captured = mean(n_captured),
  pa_detects = mean(pa_detects_total),
  tau = mean(tau),
  nsites = mean(nsites),
  ntraps = mean(ntraps),
  N_mean = mean(N_mean),
  N_mean_sd = sd(N_mean),
  N_bias = mean(N_bias),
  N_sd = mean(N_sd),
  sd_bias = (sd(N_mean) - mean(N_sd)) / mean(N_sd),
  N_naive_se = mean(N_naive_se),
  credible_proportion = sum(cred_contains_actual_N, na.rm = TRUE) / reps
), by = sim]

# Output summary statistics and all_data
write.csv(summary_stats, file = paste0("summary_stats_", directory, ".csv"))
write.csv(all_data, file = paste0("all_data_", directory, ".csv"))

toc()
