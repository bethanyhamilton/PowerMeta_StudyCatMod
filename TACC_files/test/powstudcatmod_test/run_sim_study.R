args <- commandArgs(trailingOnly = TRUE)

# Load the CSV file containing PYL_ID and the associated values
all_params <- read.csv("/home/r-environment/pyl_id_values_test.csv")  

#-------------------------------------------------------------------------------
# Parse command line argument


if (length(args) < 3) {
  stop("PYL_ID and output file path need to be arguments")
}

# Extract PYL_ID from command-line arguments
pyl_id <- as.integer(args[2])  # argument is in the form "PYL_ID #"

# Extract output file path from the third argument
output_path <- args[3]  # Argument is the file path for saving the result


#-------------------------------------------------------------------------------
# Source packages and functions


library(dplyr)
library(purrr)
library(tidyr)
library(metafor)
library(clubSandwich)
library(mvtnorm)


source("/home/r-environment/functionsdatagen.R")


#-------------------------------------------------------------------------------
# Load experimental design parameters

empirical_dat <- readRDS("/home/r-environment/dat_kjN_mathdat.rds")


#-------------------------------------------------------------------------------
# run simulations for specified batch
res <- subset(all_params, PYL_ID == pyl_id)

#res$batch <- NULL
res$PYL_ID <- NULL



tm <- system.time(res$res <- pmap(res, .f = run_sim2, 
                pilot_data = empirical_dat,
                sigma_j_sq_inc = FALSE))




#-------------------------------------------------------------------------------
# Save results and details

res$run_date <- date()
res$time <- as.numeric(tm[3])

file_name <- paste0("simulation_results_condition_test_", pyl_id, ".rds")

full_output_path <- file.path(output_path, file_name)

saveRDS(res, file = full_output_path)

