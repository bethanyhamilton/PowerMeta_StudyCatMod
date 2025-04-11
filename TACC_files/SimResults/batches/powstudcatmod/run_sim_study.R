args <- commandArgs(trailingOnly = TRUE)

#-------------------------------------------------------------------------------
# Parse command line argument


if (length(args) < 4) {
  stop("PYL_ID, output file path, and batch need to be arguments")
}

# Extract PYL_ID from command-line arguments
pyl_id <- as.integer(args[2])  # argument is in the form "PYL_ID #"

# Extract output file path from the third argument
output_path <- args[3]  # Argument is the file path for saving the result

# Extract batch number from the fourth argument
batch_file <- as.integer(args[4])  


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

# Generate list of conditions containing PYL_ID and the associated values

set.seed(03242025)

design_factors <- list(
  J = c(24, 36, 48, 60, 72),
  tau_sq = c(0.05,  0.40)^2, 
  omega_sq = c(0.05,  0.20)^2,
  rho = c(.2,  .8),
  P = c(0.05, 0.2, 0.4, 0.6, 0.8, 0.9),
  f_c_val = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"),
  bal = c("balanced_j", "unbalanced_j") 
)


batches <- 4
total_reps <- 2500

params <- expand.grid(c(design_factors, list(batch = 1:batches)))

params <- params |>
  mutate(iteration = total_reps/batches,
         seed = round(runif(1) * 2^30) + 1:n()
  ) |> 
  group_by(batch) |>
  mutate(PYL_ID = row_number()) |> 
  ungroup() |> 
  as_tibble()

params2 <- params |> filter(batch == batch_file)
params2$batch <- NULL

#-------------------------------------------------------------------------------
# run simulations for specified batch
res <- subset(params2, PYL_ID == pyl_id)

#res$batch <- NULL
res$PYL_ID <- NULL



tm <- system.time(res$res <- pmap(res, .f = run_sim, 
                pilot_data = empirical_dat,
                sigma_j_sq_inc = TRUE))

#-------------------------------------------------------------------------------
# Save results and details

res$run_date <- date()
res$time <- as.numeric(tm[3])

file_name <- paste0("simulation_results_condition_", batch_file,"_", pyl_id, ".rds")

full_output_path <- file.path(output_path, file_name)

saveRDS(res, file = full_output_path)

