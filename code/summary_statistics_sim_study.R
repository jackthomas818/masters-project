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

# directory <- "SimStudy5"

# parameter matrix from generate_parameter_matrix.R
params_matrix <- read.csv(paste0(directory, "/parameter_matrix.csv"))

# Getting model statistics
# construct summary statistics table
N_mean <- array(data = NA, dim = nrow(params_matrix))
N_lower_95 <- array(data = NA, dim = nrow(params_matrix))
N_upper_95 <- array(data = NA, dim = nrow(params_matrix))
n_captured <- array(data = NA, dim = nrow(params_matrix))
N_sd <- array(data = NA, dim = nrow(params_matrix))
N_naive_se <- array(data = NA, dim = nrow(params_matrix))

pa_detects_total <- array(data = NA, dim = nrow(params_matrix))

for (task_id in 1:nrow(params_matrix)) {
  if (task_id %% 100 == 0) {
    print(paste0(task_id, " out of ", nrow(params_matrix)))
  }

  model_file <- paste0(directory, "/model_statistics_", task_id, ".txt")
  if (file.exists(model_file)) {
    model_data <- tryCatch(read.table(model_file), error = function(e) NULL)
  }

  cr_file <- paste0(directory, "/capture-recapture-data-", task_id, ".txt")
  if (file.exists(cr_file)) {
    cr_data <- tryCatch(read.table(cr_file), error = function(e) NULL)
  }

  pa_file <- paste0(directory, "/presence-absence-data-", task_id, ".csv")
  if (file.exists(pa_file)) {
    pa_data <- tryCatch(read.csv(pa_file, header = FALSE), error = function(e) NULL)
  }

  # count number of presence detections
  pa_detects_total[task_id] <- sum(rowSums(pa_data))

  N_mean[task_id] <- unlist(model_data["Mean"])[1]
  N_lower_95[task_id] <- unlist(model_data["X2.5."])[1]
  N_upper_95[task_id] <- unlist(model_data["X97.5."])[1]
  N_sd[task_id] <- unlist(model_data["SD"])[1]
  N_naive_se[task_id] <- unlist(model_data["Naive.SE"])[1]

  if (is.null(cr_data)) {
    n_captured[task_id] <- 0
  } else {
    n_captured[task_id] <- nrow(cr_data)
  }
}

all_data <- cbind(params_matrix, N_mean, n_captured, pa_detects_total, N_lower_95, N_upper_95, N_sd, N_naive_se)

# Calculate means by simulation number
setDT(all_data)

summary_stats <- all_data[, list(
  N_actual = mean(N),
  n_captured = mean(n_captured),
  pa_detects = mean(pa_detects_total),
  tau = mean(tau),
  N_mean = mean(N_mean),
  N_sd = mean(N_sd),
  N_naive_se = mean(N_naive_se),
  lower_95_N = mean(N_lower_95),
  upper_95_N = mean(N_upper_95)
), by = sim]


# Output summary statistics and all_data
write.csv(summary_stats, file = paste0("summary_stats_", directory, ".csv"))
write.csv(all_data, file = paste0("all_data_", directory, ".csv"))

toc()
