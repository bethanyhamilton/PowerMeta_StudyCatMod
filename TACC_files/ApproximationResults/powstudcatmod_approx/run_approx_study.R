args <- commandArgs(trailingOnly = TRUE)

#-------------------------------------------------------------------------------
# Parse command line argument


if (length(args) < 3) {
  stop("PYL_ID, output, and file path need to be arguments")
}

# Extract PYL_ID from command-line arguments
pyl_id <- as.integer(args[2])  # argument is in the form "PYL_ID #"

# Extract output file path from the third argument
output_path <- args[3]  # Argument is the file path for saving the result


#-------------------------------------------------------------------------------
# Source packages and functions

source("/home/r-environment/functionsdatagen.R")

library(dplyr)
library(purrr)
library(tidyr)
library(mvtnorm)
library(stringr)




#-------------------------------------------------------------------------------
# Load experimental design parameters

empirical_dat <- readRDS("/home/r-environment/dat_kjN_mathdat.rds")

shape_rate <- MASS::fitdistr(empirical_dat$N, "gamma")
shape_rate2 <- MASS::fitdistr(empirical_dat$sigma_j_sq, "gamma")

#-------------------------------------------------------------------------------

set.seed(2503025)

design_factors <- list(
  J = c(24, 36, 48, 60, 72),
  tau_sq = c(0.05,  0.40)^2, 
  omega_sq = c(0.05,  0.20)^2,
  rho = c(.2,  .8),
  P = c(0.05, 0.2, 0.4, 0.6, 0.8, 0.9),
  f_c_val = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"),
  bal = c("balanced_j", "unbalanced_j")
)


params <- 
  expand.grid(design_factors)

params <- params |>
  dplyr::mutate(
         seed = round(runif(1) * 2^30) + 1:n()
  ) |>
  dplyr::mutate(PYL_ID = row_number())


#-------------------------------------------------------------------------------
# run simulations for specified batch
res <- subset(params, PYL_ID == pyl_id)
res$PYL_ID <- NULL


tm <- system.time(res$res <- pmap(res, .f = power_approximation,
                                  N_mean = mean(empirical_dat$N),
                                  k_mean =  mean(empirical_dat$kj),
                                  sigma_j_sq_mean = mean(empirical_dat$sigma_j_sq),
                                  N_dist = shape_rate, 
                                  sigma_j_sq_dist = shape_rate2,
                                  pilot_data = empirical_dat, 
                                  average_power = TRUE,
                                  iterations = 100,
                                  sample_size_method = c("balanced","stylized","empirical"))) 
  


#-------------------------------------------------------------------------------
# Save results and details

res$run_date <- date()
res$time <- as.numeric(tm[3])

file_name <- paste0("approx_results_condition_", pyl_id, ".rds")

full_output_path <- file.path(output_path, file_name)

saveRDS(res, file = full_output_path)







