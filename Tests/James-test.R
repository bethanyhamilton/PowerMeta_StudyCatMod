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
# Test generate_meta()

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

# Generate big dataset with ES parameters

multi <- 20
sample_dat <- n_ES_empirical(empirical_dat, J = multi * J)

big_param_dat <- generate_meta(
  J = multi * J, 
  tau_sq = tau_sq, 
  omega_sq = omega_sq, 
  bal = bal, 
  rho = rho, 
  sample_sizes = sample_dat$N,
  k_j = sample_dat$kj,
  mu_vector = mu_vector,
  f_c_val = f_c_val,
  return_study_params = TRUE
) %>% 
  unnest(c(delta, X))

# Compare generating parameters to estimates

library(lme4)
ES_param_fit <- lmer(delta ~ 0 + category + (1 | studyid), data = big_param_dat)
summary(ES_param_fit)
data.frame(param = mu_vector, est = fixef(ES_param_fit))
data.frame(param = sqrt(c(tau_sq, omega_sq)), est = VarCorr(ES_param_fit))

# Generate big dataset with ES estimates

big_dat <- generate_meta(
  J = multi * J, 
  tau_sq = tau_sq, 
  omega_sq = omega_sq, 
  bal = bal, 
  rho = rho, 
  sample_sizes = sample_dat$N, 
  k_j = sample_dat$kj,
  mu_vector = mu_vector,
  f_c_val = f_c_val,
  return_study_params = FALSE
) %>%
  group_by(studyid) %>%
  mutate(var_g_j = mean(var_g, na.rm = TRUE))
  

# Compare generating parameters to estimates

library(metafor)

V_list <- 
  vcalc(
    vi = var_g_j,
    cluster = studyid,
    rho = rho,
    obs = esid,
    sparse = TRUE,
    data = big_dat 
  )

rma_fit <- rma.mv(
  g ~ 0 + category,
  V = V_list,
  random = ~ 1 | studyid / esid,
  data = big_dat,
  sparse = TRUE
)
summary(rma_fit)
data.frame(param = mu_vector, est = coef(rma_fit))
data.frame(param = c(tau_sq, omega_sq), est = rma_fit$sigma2) |> sqrt()


#-------------------------------------------------------------------------------
# Test run_power()
J <- 48
run_power(
  C = 3, 
  J = J, 
  tau_sq = .05^2, 
  omega_sq = .05^2, 
  rho = 0.6, 
  k_j = 1L + rpois(J, 3),
  P = 10000, 
  f_c_val = "P4",
  mu_vector = c(0,0.1,0.2),
  bal = "balanced_j",
  sigma_j_sq = rep(4 / 60, J)
)

run_power(
  C = 3, 
  J = J, 
  tau_sq = .05^2, 
  omega_sq = .05^2, 
  rho = 0.6, 
  k_j = rep(1L, J),
  P = 10000, 
  f_c_val = "P4",
  mu_vector = c(0,0.1,0.2),
  bal = "balanced_j",
  sigma_j_sq = rep(4 / 60, J)
)


debug(run_power)
run_power(
  C = length(mu_vector), 
  J = J, 
  tau_sq = tau_sq, 
  omega_sq = omega_sq, 
  rho = rho, 
  k_j = param_dat$k_j,
  P = 0.4, 
  f_c_val = f_c_val,
  mu_vector = mu_vector,
  bal = bal,
  sigma_j_sq = param_dat$var_g_j
)

# generate study-level features
sample_dat <- n_ES_empirical(empirical_dat, J = J)

param_dat <- generate_meta(
  J = J, 
  tau_sq = tau_sq, 
  omega_sq = omega_sq, 
  bal = bal, 
  rho = rho, 
  sample_sizes = sample_dat$N,
  k_j = sample_dat$kj,
  mu_vector = mu_vector,
  f_c_val = f_c_val,
  return_study_params = TRUE,
  seed = 20250430
) %>% 
  unnest(c(delta, X)) %>%
  group_by(studyid) %>%
  summarize(across(c(k_j,N,category), unique))


# Generate dataset with ES estimates

ES_dat <- generate_meta(
  J = J, 
  tau_sq = tau_sq, 
  omega_sq = omega_sq, 
  bal = bal, 
  rho = rho, 
  sample_sizes = sample_dat$N, 
  k_j = sample_dat$kj,
  mu_vector = mu_vector,
  f_c_val = f_c_val,
  seed = 20250430
) %>%
  group_by(studyid) %>%
  mutate(var_g_j = mean(var_g, na.rm = TRUE))

param_dat <- 
  ES_dat %>%
  group_by(studyid) %>%
  summarize(var_g_j = unique(var_g_j)) %>%
  inner_join(param_dat, by = "studyid")

# Fit CHE with known variance components
V_list <- 
  vcalc(
    vi = var_g_j,
    cluster = studyid,
    rho = rho,
    obs = esid,
    sparse = TRUE,
    data = ES_dat 
  )

rma_fit <- rma.mv(
  g ~ 0 + category,
  V = V_list,
  random = ~ 1 | studyid / esid,
  data = ES_dat,
  sigma2 = c(tau_sq, omega_sq),
  sparse = TRUE
)

res_RVE <- conf_int(rma_fit, vcov = "CR2")
SEs <- res_RVE$SE
df <- res_RVE$df
Wald_stat <- Wald_test(rma_fit, constraints = constrain_equal(1:4), vcov = "CR2")
Wald_stat

cat_res <- multiple_categories(
  data = param_dat, 
  moderator_val = param_dat$category, 
  cluster_id = param_dat$studyid, 
  sigma_j_sq_val = param_dat$var_g_j, 
  rho_val = rho, 
  omega_sq_val = omega_sq, 
  tau_sq_val = tau_sq
)
cat_res
res_RVE
cat_res$W
1 / res_RVE$SE^2

power_CHE_RVE_study_cat(
  data = param_dat, 
  moderator_val = param_dat$category, 
  cluster_id = param_dat$studyid, 
  sigma_j_sq_val = 4 / param_dat$N, 
  rho_val = rho, 
  omega_sq_val = omega_sq, 
  tau_sq_val = tau_sq,
  mu = mu_vector, 
  alpha = .05
)
