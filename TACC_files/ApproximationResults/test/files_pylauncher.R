#write pylauncher command lines files (https://github.com/TACC/pylauncher,
#https://docs.tacc.utexas.edu/software/pylauncher/)

library(dplyr)

set.seed(04222025)


#8640 conditions
design_factors <- list(
  J = c(24, 36, 48, 60, 72),
  tau_sq = c(0.05, 0.20, 0.40)^2, 
  omega_sq = c(0.05, 0.20)^2,
  rho = c(.2, 0.5, .8),
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

# run only the conditions with tau_sq = 0.0025
params <- params |> filter(tau_sq == "0.0025")


# Open the file for writing (no extension)
FileName2 <- paste("TACC_files/ApproximationResults/test/command_lines_approx_test", sep="")

# Open the file for writing
file_conn <- file(FileName2, "w")


for (i in 1:dim(params)[1]) {
  # Construct the command string
  command <- paste('apptainer run --bind /scratch/08147/bethanyh/pylauncher_out_approx_test:/home/r-environment/output /scratch/08147/bethanyh/powstudcatmod_approx_v0_test.sif R -e "source(\'/home/r-environment/run_approx_study.R\')" --args PYL_ID ', params$PYL_ID[i], ' "/home/r-environment/output" ' ,sep = "")

  # Write the command to the file
  writeLines(command, file_conn)
}


# Close the connection to the file
close(file_conn)

