#-------------------------------------------------------------------------------
# Source packages and functions


library(dplyr)
library(purrr)
library(tidyr)
library(metafor)
library(clubSandwich)
library(mvtnorm)
library(stringr)

source("SimFunctions/functionsdatagen.R")
#-------------------------------------------------------------------------------
# Load experimental design parameters

# empirical data sets clustered at study level
empirical_dat <- readRDS("SimFunctions/dat_kjN_mathdat.rds")
shape_rate <- MASS::fitdistr(empirical_dat$N, "gamma")
shape_rate2 <- MASS::fitdistr(empirical_dat$sigma_j_sq, "gamma")


# empirical data clustered at sample level
empirical_dat_sample <- readRDS("SimFunctions/dat_kjN_mathdat_alt.rds")
shape_rate_sample <- MASS::fitdistr(empirical_dat_sample$N, "gamma")
shape_rate2_sample <- MASS::fitdistr(empirical_dat_sample$sigma_j_sq, "gamma")

#-------------------------------------------------------------------------------
set.seed(20250429)


design_factors <- list(
  J = c(48),
  tau_sq = c(0.05)^2, 
  omega_sq = c(0.20)^2,
  rho = c(.2),
  P = c(0.4),
  f_c_val = c("P8"),
  bal = c("balanced_j") 
)


batches <- 1
total_reps <- 100

params <- expand.grid(c(design_factors, list(batch = 1:batches)))

params <- params |>
  mutate(iteration = total_reps/batches,
         seed = round(runif(1) * 2^30) + 1:n()
  ) |> 
  group_by(batch) |>
  ungroup() |> 
  as_tibble()


#-------------------------------------------------------------------------------
res <- params

res$batch <- NULL

res_samplelvl <- res


tm <- system.time(res$res <- pmap(res, .f = run_sim_alt, 
                pilot_data = empirical_dat,
                N_dist = shape_rate, 
                sigma_j_sq_dist = shape_rate2,
                sample_size_method = c("balanced","stylized","empirical"),
                sigma_j_sq_inc = TRUE))


tm2 <- system.time(res_samplelvl$res <- pmap(res_samplelvl, .f = run_sim_alt, 
                                  pilot_data = empirical_dat_sample,
                                  N_dist = shape_rate_sample, 
                                  sigma_j_sq_dist = shape_rate2_sample,
                                  sample_size_method = c("balanced","stylized","empirical"),
                                  sigma_j_sq_inc = TRUE))


#-------------------------------------------------------------------------------
# Save results and details

# study lvl empirical dat results
res$run_date <- date()
res$time <- as.numeric(tm[3])

saveRDS(res, file = "Tests/simulation_results_condition_test_study.rds")



# sample lvl empirical dat results
res_samplelvl$run_date <- date()
res_samplelvl$time <- as.numeric(tm2[3])

saveRDS(res_samplelvl, file = "Tests/simulation_results_condition_test_sample.rds")


#-------------------------------------------------------------------------------
# Analyze Results

 results_studylvl <- res |> 
  unnest(res, names_sep = ".") |> 
  dplyr::summarise(
    n_iterations_conv = sum(!is.na(res.df1)),
    n_iterations_noconv = sum(is.na(res.df1)),
    n_iterations = n(),
    conv_rate =  mean(n_iterations_conv/100),
    p_val_avg = mean(res.p_val),
    df2_mean = mean(res.df2),
    rej_rate_05 = mean(res.p_val < 0.05, na.rm = TRUE),
    MCSE_rej_rate_05 = sqrt(rej_rate_05*(1 - rej_rate_05)/n_iterations_conv),
    mu_est_mean = list(reduce(res.est, `+`) / n()),
    mu_est_var = list(map(transpose(res.est), ~ var(unlist(.)))),
    mu_params = list(reduce(res.mu_vector_list, `+`) / n()),
    mu_est_bias = map2(mu_est_mean, mu_params, ~ .x - .y), 
    mean_power_empirical = mean(res.power_empirical),
    se_power_empirical = sd(res.power_empirical)/sqrt(n_iterations),
    power_app_sim_diff =  mean_power_empirical - rej_rate_05,
    .groups = "drop"
  )  

results_studylvl$power_app_sim_diff


results_samplelvl <- res_samplelvl |> 
  unnest(res, names_sep = ".") |> 
  dplyr::summarise(
    n_iterations_conv = sum(!is.na(res.df1)),
    n_iterations_noconv = sum(is.na(res.df1)),
    n_iterations = n(),
    conv_rate =  mean(n_iterations_conv/100),
    p_val_avg = mean(res.p_val),
    df2_mean = mean(res.df2),
    rej_rate_05 = mean(res.p_val < 0.05, na.rm = TRUE),
    MCSE_rej_rate_05 = sqrt(rej_rate_05*(1 - rej_rate_05)/n_iterations_conv),
    mu_est_mean = list(reduce(res.est, `+`) / n()),
    mu_est_var = list(map(transpose(res.est), ~ var(unlist(.)))),
    mu_params = list(reduce(res.mu_vector_list, `+`) / n()),
    mu_est_bias = map2(mu_est_mean, mu_params, ~ .x - .y), 
    mean_power_empirical = mean(res.power_empirical),
    se_power_empirical = sd(res.power_empirical)/sqrt(n_iterations),
    power_app_sim_diff =  mean_power_empirical - rej_rate_05,
    .groups = "drop"
  )  

results_samplelvl$power_app_sim_diff
