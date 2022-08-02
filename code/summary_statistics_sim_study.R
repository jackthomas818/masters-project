# Load packages
library(readr)
library(tictoc)
library(data.table)

tic()

print("finished loading libraries")

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# directory of simulation study
directory <- args[1]

#directory <- "SimStudy4/"

# Getting presence-absence data
pa_files <- list.files(directory, pattern = "presence-absence")
pa_data <- lapply(paste0(directory, "\\", pa_files), function(i) {
  read.csv(i, header = FALSE)
})

# Getting capture-recapture data
cr_files <- list.files(directory, pattern = "capture-recapture")
cr_data <- lapply(paste0(directory, "\\", cr_files), function(i) {
  read.table(i, header = FALSE)
})

# parameter matrix from generate_parameter_matrix.R
params_matrix <- read.csv("./parameter_matrix.csv")

# Getting model statistics
# construct summary statistics table
N_mean <- array(data = NA, dim = nrow(params_matrix))
N_lower_95 <- array(data = NA, dim = nrow(params_matrix))
N_upper_95 <- array(data = NA, dim = nrow(params_matrix))
n_captured <- array(data = NA, dim = nrow(params_matrix))
N_sd <- array(data = NA, dim = nrow(params_matrix))
N_naive_se <- array(data = NA, dim = nrow(params_matrix))

for (task_id in 1:nrow(params_matrix)) {
  model_data <- read.table(paste0(directory, "/model_statistics_", task_id, ".txt"))
  print(nrow(model_data))
  cr_data <- read.table(paste0(directory, "/capture-recapture-data-", task_id, ".txt"))

  N_mean[task_id] <- unlist(model_data["Mean"])[1]
  N_lower_95[task_id] <- unlist(model_data["X2.5."])[1]
  N_upper_95[task_id] <- unlist(model_data["X97.5."])[1]
  N_sd[task_id] <- unlist(model_data["SD"])[1]
  N_naive_se[task_id] <- unlist(model_data["Naive.SE"])[1]
  n_captured[task_id] <- nrow(cr_data)
}


all_data <- cbind(params_matrix, N_mean, n_captured, N_lower_95, N_upper_95, N_sd, N_naive_se)


# Calculate means by simulation number
setDT(all_data)

all_data[, list(
  N_actual = mean(N),
  n_captured = mean(n_captured),
  tau = mean(tau),
  N_mean = mean(N_mean),
  N_sd = mean(N_sd),
  N_naive_se=mean(N_naive_se),
  lower_95_N = mean(N_lower_95),
  upper_95_N = mean(N_upper_95)
), by = sim]

toc()
