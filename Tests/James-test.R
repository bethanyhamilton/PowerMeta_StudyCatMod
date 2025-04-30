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

J <- 48
tau_sq <- 0.05^2
omega_sq = 0.20^2
rho <- 0.2
P <- 0.4
f_c_val <- "P8"
bal <- "balanced_j"

N_mean <- mean(empirical_dat$N) 
k_mean <- mean(empirical_dat$kj)
sigma_j_sq_mean <- mean(empirical_dat$sigma_j_sq)

mu_vector <- mu_values(
  J = J, tau_sq = tau_sq, omega_sq = omega_sq, 
  rho = rho, P = P,
  f_c_val = f_c_val,
  bal = bal, 
  k_j  = k_mean,  
  sigma_j_sq = sigma_j_sq_mean, 
  N = N_mean
)


sample_dat <- n_ES_empirical(empirical_dat, J = J)

debug(generate_meta)
multi <- 100

big_param_dat <- generate_meta(
  J = multi * J, 
  tau_sq = tau_sq, 
  omega_sq = omega_sq, 
  bal = bal, 
  rho = rho, 
  sample_sizes = rep(sample_dat$N, multi),
  k_j = rep(sample_dat$kj,multi),
  mu_vector = mu_vector,
  f_c_val = f_c_val,
  return_study_params = TRUE
) %>% 
  unnest(c(delta, X))

library(lme4)
lmer(delta ~ 0 + category + (1 | studyid), data = big_param_dat) |>
  summary()


big_dat <- generate_meta(
  J = multi * J, 
  tau_sq = tau_sq, 
  omega_sq = omega_sq, 
  bal = bal, 
  rho = rho, 
  sample_sizes = rep(sample_dat$N, multi),
  k_j = rep(sample_dat$kj,multi),
  mu_vector = mu_vector,
  f_c_val = f_c_val,
  return_study_params = FALSE
)

library(metafor)
rma.mv()