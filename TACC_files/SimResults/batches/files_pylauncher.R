#write pylauncher command lines files (https://github.com/TACC/pylauncher,
#https://docs.tacc.utexas.edu/software/pylauncher/)

library(dplyr)

set.seed(20250422)

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

batch_file <-  4
params2 <- params %>% filter(batch == batch_file)
params2$batch <- NULL



# Open the file for writing (no extension)
FileName2 <- paste("TACC_files/SimResults/batches/command_lines",batch_file,sep="")



# Open the file for writing
file_conn <- file(FileName2, "w")

# Loop through the numbers 1 to 3840 and write the commands to the file
# for (i in 1:dim(params2)[1]) {
#   # Construct the command string
#   command <- paste("Rscript run_sim_study.R PYL_ID", i)
# 
#   # Write the command to the file
#   writeLines(command, file_conn)
# }



for (i in 1:dim(params2)[1]) {
  # Construct the command string
  command <- paste('apptainer run --bind /scratch/08147/bethanyh/pylauncher_out', batch_file,':/home/r-environment/output /scratch/08147/bethanyh/powstudcatmod_v0.sif R -e "source(\'/home/r-environment/run_sim_study.R\')" --args PYL_ID ', params2$PYL_ID[i], ' "/home/r-environment/output" ', batch_file ,sep = "")

  # Write the command to the file
  writeLines(command, file_conn)
}


# Close the connection to the file
close(file_conn)

