rm(list=ls())
library(tidyverse)

# processing simulation and approximation results

#-------------------------------------------------------------------------------
## read in simulation results -- raw data. Four batches each with 
## one file for each condition (of 635 iterations) due to parametric job launcher.
## across four batches the condition should have 2,500 iterations.
### takes a bit to run...

## to reduce bottleneck --  read in simulation results and summarize
## get all names of rds files > nest by condition > read in the four files for 
## each condition > and then summarize

setwd("dat/sim/pylauncher_out")


#function to read in and summarize sim files. due to non-convergence, each condition 
#may not have the same number of samples, so pass a value to slice argument 
#to make sure that each condition has the same number of samples.
process_sim_files <- function(files, summary = TRUE, slice = NULL){
  
  sim_results_raw <- files |>
    map(readRDS) |>
    bind_rows() |>
    unnest(res, names_sep = ".")
  
  nonconverged <- sim_results_raw |>
    filter(is.na(res.df1)) 
  
  
  if(summary){
    
    if(!is.null(slice)) {
      summary_converged <-  sim_results_raw |> 
        filter(!is.na(res.df1)) |> 
        mutate(n_iterations_conv = n()) |> 
        slice(1:slice) |> 
        group_by(J, tau_sq, omega_sq, rho, P, f_c_val, bal) |> 
        dplyr::summarise(
          n_iterations_retained = n(),
          conv_rate =  mean(n_iterations_conv/2500),
          p_val_avg = mean(res.p_val),
          df2_mean = mean(res.df2),
          rej_rate_05 = mean(res.p_val < 0.05, na.rm = TRUE),
          MCSE_rej_rate_05 = sqrt(rej_rate_05*(1 - rej_rate_05)/n_iterations_retained),
          mu_est_mean = list(reduce(res.est, `+`) / n()),
          mu_est_var = list(map(transpose(res.est), ~ var(unlist(.)))),
          mu_params = list(reduce(res.mu_vector_list, `+`) / n()),
          mu_est_bias = map2(mu_est_mean, mu_params, ~ .x - .y), 
          .groups = "drop"
        ) 
    }else{
      summary_converged <-  sim_results_raw |> 
        filter(!is.na(res.df1)) |> 
        group_by(J, tau_sq, omega_sq, rho, P, f_c_val, bal) |> 
        dplyr::summarise(
          n_iterations_conv = n(),
          conv_rate = mean(n_iterations_conv/2500),
          p_val_avg = mean(res.p_val),
          df2_mean = mean(res.df2),
          rej_rate_05 = mean(res.p_val < 0.05, na.rm = TRUE),
          MCSE_rej_rate_05 = sqrt(rej_rate_05*(1 - rej_rate_05)/n_iterations_conv),
          mu_est_mean = list(reduce(res.est, `+`) / n()),
          mu_est_var = list(map(transpose(res.est), ~ var(unlist(.)))),
          mu_params = list(reduce(res.mu_vector_list, `+`) / n()),
          mu_est_bias = map2(mu_est_mean, mu_params, ~ .x - .y), 
          
          .groups = "drop"
        ) 
    }
    ## maybe add df of the sim results later...example below from previous study
    # est_V = var(est, na.rm = TRUE), #done
    # var_M = mean(est_var, na.rm = TRUE),
    # var_df = 2 * var_M^2 / var(est_var, na.rm = TRUE)
    
    return(list(converged_results = summary_converged, nonconverged_results = nonconverged))
    
  } else{
    
    converged <-  sim_results_raw |> 
      filter(!is.na(res.df1))
    
    return(list(converged_results = converged, nonconverged_results = nonconverged))
  }

}


sim_data  <- list.files(pattern = ".rds") |> 
  tibble() |> 
  rename(files = `list.files(pattern = ".rds")`) |> 
  mutate(condition_number = str_extract(files, "\\d+(?=\\.rds$)")) |> 
  nest(data= -condition_number) |> 
  mutate(
    processed  = map(data, ~ map(.x, process_sim_files, summary = TRUE, slice = 2496)),
    #processed  = map(data, ~ map(.x, process_sim_files, summary = TRUE)),
    converged = map(processed, ~ map_dfr(.x, "converged_results")),
    nonconverged = map(processed, ~ map_dfr(.x, "nonconverged_results"))) |> 
  select(condition_number, converged, nonconverged)

sim_results <- sim_data |> 
  select(condition_number, converged) |>  
  unnest(converged)


sim_results_nonconvergence <- sim_data |> 
  select(condition_number, nonconverged) |>  
  unnest(nonconverged)


save(sim_results, file = "../../../Analysis/sim_results.RData")

save(sim_results_nonconvergence, file = "../../../Analysis/sim_resultsnonconvergence.RData")

#-------------------------------------------------------------------------------
## approximation results

setwd("../../approx/pylauncher_out_approx/")

approx_results <- list.files(pattern = ".rds") |> 
  map(readRDS) |>  
  bind_rows() |> 
  unnest(res, names_sep = ".") |> 
  pivot_wider(
    names_from = res.samp_method,
    values_from = c(res.power.mean, res.ncp.mean, res.df_den.mean, 
                    res.power.var, res.ncp.var, res.df_den.var, 
                    res.power.n, res.ncp.n, res.df_den.n,
                    res.power.se, res.ncp.se, res.df_den.se)
  ) 




save(approx_results, file = "../../../Analysis/approx_results.RData")


