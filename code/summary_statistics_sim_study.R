# Load packages
library(readr)
library(tictoc)

tic()

print("finished loading libraries")

# get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# directory of simulation study
directory <- args[1]

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

# Getting model statistics
model_files <- list.files(directory, pattern = "model_statistics")
model_data <- lapply(paste0(directory, "\\", model_files), function(i) {
  read.table(i, header = TRUE)
})

toc()
