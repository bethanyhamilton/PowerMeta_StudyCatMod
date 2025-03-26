rm(list=ls())
source("SimFunctions/functionsdatagen.R")

library(dplyr)
library(purrr)
library(tidyr)
library(mvtnorm)
library(plyr)
library(future)
library(furrr)


plan(multisession)
set.seed(2503025)

dat_kjN <- readRDS("SimFunctions/dat_kjN_mathdat.rds")

shape_rate <- MASS::fitdistr(dat_kjN$N, "gamma")
shape_rate2 <- MASS::fitdistr(dat_kjN$sigma_j_sq, "gamma")

design_factors <- list(
  J = c(24, 36, 48, 60, 72),
  tau_sq = c(0.05,  0.40)^2, 
  omega_sq = c(0.05,  0.20)^2,
  rho = c(.2,  .8),
  P = c(0.05, 0.2, 0.4, 0.6, 0.8, 0.9),
  f_c_val = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"),
  bal = c("balanced_j", "unbalanced_j")
)


pow_params <- 
  expand.grid(design_factors)

pow_params <- pow_params |>
  dplyr::mutate(
         seed = round(runif(1) * 2^30) + 1:n()
  ) 



#parallel::detectCores()
#source_obj <- ls()
#cluster=snow::makeCluster(7)
#doSNOW::registerDoSNOW(cluster)

#snow::clusterExport(cluster, source_obj)
#snow::clusterApply(cluster, library("tidyverse"))

# pow_params <- pow_params |> 
#   slice(1:5)



tm <- system.time(app_results <-
  pow_params %>%
  mutate(res = future_pmap(., .f = purrr::possibly(power_approximation, otherwise = data.frame()),
                           N_mean = mean(dat_kjN$N),
                           k_mean =  mean(dat_kjN$kj),
                           sigma_j_sq_mean = mean(dat_kjN$sigma_j_sq),
                           N_dist = shape_rate, 
                           sigma_j_sq_dist = shape_rate2,
                           pilot_data = dat_kjN, 
                           average_power = TRUE,
                           iterations = 1,
                           sample_size_method = c("balanced","stylized","empirical"))) %>%
  tidyr::unnest(cols = res) )

tm
app_results

pow_params %>%
  anti_join(app_results) %>%
  nrow() %>% 
  identical(0L)


session_info <- sessionInfo()
run_date <- date()



save(tm,   pow_params, app_results, session_info, run_date, file = "Approximations/approximation_resutls.RData")



# tm <- system.time(results <- plyr::mdply(params, 
#                                          .fun = power_approximation,
#                                          N_mean = mean(dat_kjN$N),
#                                          k_mean =  mean(dat_kjN$kj),
#                                          sigma_j_sq_mean = mean(dat_kjN$sigma_j_sq),
#                                          N_dist = shape_rate, 
#                                          sigma_j_sq_dist = shape_rate2,
#                                          pilot_data = dat_kjN, 
#                                          average_power = TRUE,
#                                          iterations = 1,
#                                          sample_size_method = c("balanced","stylized","empirical"),
#                                          .inform = TRUE,
#                                          .parallel = FALSE))
# 
# tm
#parallel::stopCluster(cluster)






