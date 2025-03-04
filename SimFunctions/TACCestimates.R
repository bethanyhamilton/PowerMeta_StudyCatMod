library(tidyverse)
library(metafor)
library(clubSandwich)
library(mvtnorm)

rm(list=ls())

source("SimFunctions/functionsdatagen.R")

design_factors <- list(
  J = c(24, 36, 48, 60, 72),
  tau_sq = c(0.05,  0.40)^2, 
  omega_sq = c(0.05,  0.20)^2,
  rho = c(.2,  .8),
  P = c(0.05, 0.2, 0.4, 0.6, 0.8, 0.9),
  f_c_val = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"),
  bal = c("balanced_j", "unbalanced_j")
)


params <- expand.grid(design_factors) %>% 
  mutate(
    iterations = 20,
    seed = round(runif(1) * 2^30) + 1:n()
  ) %>%
  as_tibble()

set.seed(0303202571)
 params <- params |> 
   group_by(J) |> 
   sample_n(size=10, replace = TRUE) 



 empirical_dat <- readRDS("SimFunctions/dat_kjN_mathdat.rds")
 
 

tm <- system.time(  results <- plyr::mdply(params,
                         .fun = run_sim2,
                         pilot_data = empirical_dat,
                         sigma_j_sq_inc = FALSE,
                         .parallel = FALSE))

sum(results$time)
tm


#saveRDS(results, "SimFunctions/TACCtest.rds")

dat_kjN <- readRDS("SimFunctions/TACCtest.rds")

tm_j<- results |> 
  group_by(J) |> 
  summarize(mean_tm_20rep = mean(time),
            avg_tm_per_cond = mean_tm_20rep*(2500/20),
            tot_tm_byJ = avg_tm_per_cond*768) 

tm_j 

# total compute time in hrs
as.numeric(sum(tm_j$tot_tm_byJ))/60^2

  
  # params2 <- params %>% filter(seed == 245717973)
  # results <- plyr::mdply(params2,
  #                        .fun = run_sim2,
  #                        pilot_data = empirical_dat,
  #                        sigma_j_sq_inc = FALSE,
  #                        .parallel = FALSE)
  # 
  # set.seed(245717973)
  # sample_dat <- n_ES_empirical(empirical_dat, J = 36)
  
