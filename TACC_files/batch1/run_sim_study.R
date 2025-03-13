args <- commandArgs(trailingOnly = TRUE)

# Load the CSV file containing PYL_ID and the associated values
all_params <- read.csv("TACC_files/pyl_id_values1.csv")  

#-------------------------------------------------------------------------------
# Parse command line argument

# row_to_run <- 
#   args |>
#   paste(collapse = " ") |>
#   stringr::str_extract("batch [0-9]+") |>
#   stringr::str_sub(7, -1) |>
#   as.integer()


if (length(args) < 2) {
  stop("PYL_ID needs to be an argument")
}

# Extract PYL_ID from command-line arguments
pyl_id <- as.integer(args[2])  # argument is in the form "PYL_ID #"


#-------------------------------------------------------------------------------
# Source packages and functions


# require(devtools)
# install_version("metafor", version = "4.8.0", repos = "http://cran.us.r-project.org")
# library(metafor)

# if(!require(pacman)){
#   install.packages("pacman")
# }
# 
# p_load(
#   tidyverse, 
#   clubSandwich,
#   metafor,
#   mvtnorm,
#   purrr,
#   future,
#   furrr#,
# #  snow,
# #  doSNOW
# )


library(tidyverse)
library(metafor)
library(clubSandwich)
library(mvtnorm)


source("SimFunctions/functionsdatagen.R")
#source("SimFunctions/run_sim.R")

#-------------------------------------------------------------------------------
# Load experimental design parameters

empirical_dat <- readRDS("SimFunctions/dat_kjN_mathdat.rds")


#-------------------------------------------------------------------------------
# run simulations for specified batch
res <- subset(all_params, PYL_ID == pyl_id)

#res$batch <- NULL
res$PYL_ID <- NULL

tic()

res$res <- pmap(res, .f = run_sim2, 
                pilot_data = empirical_dat,
                sigma_j_sq_inc = FALSE)

tm <- toc(quiet = TRUE)


#-------------------------------------------------------------------------------
# Save results and details

res$run_date <- date()
res$time <- tm$toc - tm$tic

# maybe should add which batch to save file name as well. 
saveRDS(res, file = paste0("simulation_results_condition", pyl_id, ".rds"))