library(dplyr)
library(purrr)
library(tidyr)
library(metafor)
library(clubSandwich)
library(mvtnorm)
library(stringr)


#######Parameters 
# J = number of studies
# tau_sq = between-study variance
# omega_sq = within-study variance
# rho = sample correlation
# P = power level
# k_j = number of effect sizes per study
# f_c_val = pattern of the mu values 
# bal = the balance of the number of studies across categories (currently only balanced)
# sigma_j_sq = the average sample variance of effect sizes at each study (can either specify sigma_j_sq or N)
# N = primary study sample size (can either specify sigma_j_sq or N)


#------------------------------------------------------------------------------------
### ------------------------------------------- ###
###         Test Power and mu vector            ###
### ------------------------------------------- ###

## GOAL: see if I get same power as I used to generate mu when I plug in mu values
## was able to get 0.9 again

rm(list=ls())

source("SimFunctions/functionsdatagen.R")

## make my mu vector
mu_vector <- mu_values(J = 12, tau_sq = .05^2, 
                       omega_sq = .05^2,
                       rho = .5, P = .9, 
                       k_j = 3, N = 30, 
                       f_c_val = "P5", 
                       #sigma_j_sq = NA,
                       bal ="balanced_j" )
mu_vector


run_power(C = 4,
          J = 12,
          tau_sq = .05^2,
          omega_sq = .05^2,
          rho = .5,
          k_j = 3,
          N = 30,
          bal = "balanced_j",
          P = .9, 
          f_c_val = "P5",
          mu_vector = mu_vector)



### --------------------------------------------------------- ###
### Demonstrate Power Approximation given a set of conditions ###
### --------------------------------------------------------- ###
rm(list=ls())
source("SimFunctions/functionsdatagen.R")
set.seed(2202025)

dat_kjN <- readRDS("SimFunctions/dat_kjN_mathdat.rds")
dat_kjN_samp <- n_ES_empirical(dat_kjN, J = 24, with_replacement = TRUE)
shape_rate <- MASS::fitdistr(dat_kjN$N, "gamma")
shape_rate2 <- MASS::fitdistr(dat_kjN$sigma_j_sq, "gamma")


test <- power_approximation(
  J = 24, 
  tau_sq = .40^2, 
  omega_sq= .10^2, 
  rho = 0.8, 
  N_mean = mean(dat_kjN$N), 
  k_mean = mean(dat_kjN$kj),
  sigma_j_sq_mean = mean(dat_kjN$sigma_j_sq),
  sigma_j_sq_dist = shape_rate2,
  N_dist = shape_rate, 
  pilot_data = dat_kjN_samp, 
  iterations = 10,
  sample_size_method = c("balanced","stylized","empirical"),
  P = .5,
  f_c_val = "P5",
  bal = "balanced_j",
  average_power = FALSE,
    
  seed = NULL)
test

tm <- system.time(blah <- power_approximation(
  J = 24, 
  tau_sq = .40^2, 
  omega_sq= .10^2, 
  rho = 0.8, 
  N_mean = mean(dat_kjN$N), 
  k_mean = mean(dat_kjN$kj),
  sigma_j_sq_mean = mean(dat_kjN$sigma_j_sq),
  sigma_j_sq_dist = shape_rate2,
  N_dist = shape_rate, 
  pilot_data = dat_kjN_samp, 
  iterations = 100,
  sample_size_method = c("balanced","stylized","empirical"),
  P = .5,
  f_c_val = "P5",
  bal = "balanced_j",
  average_power = TRUE,
  seed = NULL))


#------------------------------------------------------------------------------------

### ------------------------------------------- ###
###         Test generate Meta and Df           ###
### ------------------------------------------- ###
rm(list = ls())
source("SimFunctions/functionsdatagen.R")


sample_empirical_dat <- tibble(N = rep(200, 12),
                               k_j  = c(3, 4, 3, 5, 3, 3, 3, 3, 3, 3, 3, 3))
mu_vector <- mu_values(
  J = 12,
  tau_sq = .05^2,
  omega_sq = .05^2,
  rho = .5,
  P = .9,
  k_j = sample_empirical_dat$k_j,
  N = sample_empirical_dat$N,
  f_c_val = "P5",
  #sigma_j_sq = NA,
  bal = "balanced_j"
)

meta_dat <- generate_meta(
  J = 12,
  tau_sq = .05^2,
  omega_sq = .05^2,
  bal = "balanced_j",
  mu_vector = mu_vector,
  #C = 4,
  rho = .5,
  sample_sizes = sample_empirical_dat$N,
  k_j = sample_empirical_dat$k_j,
  
  f_c_val = "P5",
  return_study_params = FALSE,
  seed = NULL
)




### clubSandwich df
V <- with(
  meta_dat,
  clubSandwich::impute_covariance_matrix(
    vi = var_g,
    cluster = studyid,
    r = .5,
    smooth_vi = TRUE
  )
)

fit <- rma.mv(
  g ~ -1 + category,
  V = V,
  random = ~ 1 | studyid / esid,
  sparse = TRUE,
  data = meta_dat,
  test = "t",
  method = "REML"
)



V_sep <- vcovCR(fit, cluster = meta_dat$studyid, type = "CR2")

clubsandwichdf <- Wald_test((fit), constraints = constrain_equal(1:4), vcov = V_sep)
#clubsandwichdf$df_denom

### my function df
multiplecat  <- multiple_categories(
  data = meta_dat,
  moderator_val = category,
  cluster_id = studyid,
  sigma_j_sq_val = var_g,
  rho_val = .5,
  omega_sq_val = fit$sigma2[2],
  tau_sq_val = fit$sigma2[1]
)

df_Z <- df_val(E_Vr = multiplecat$E_Vr,
               nu_Vr = multiplecat$nu_Vr,
               W = multiplecat$W)[[2]]
q <- df_val(E_Vr = multiplecat$E_Vr,
            nu_Vr = multiplecat$nu_Vr,
            W = multiplecat$W)[[1]]

C_sep <- constrain_equal(1:4, coefs = coef(fit))
mu_hat <- coef(fit)
VR <- vcovCR(fit, type = "CR2") |> as.matrix()


Q <- as.numeric(t(C_sep %*% mu_hat) %*% solve(C_sep %*% VR %*% t(C_sep)) %*% (C_sep %*% mu_hat))

delta <- (df_Z - q + 1) / (df_Z)
Fstat <- delta * Q / q
df_num <- q
df_den <- df_Z - q + 1


## test
all.equal(df_num, clubsandwichdf$df_num)
all.equal(df_den, clubsandwichdf$df_denom)
all.equal(delta, clubsandwichdf$delta)
all.equal(Fstat, clubsandwichdf$Fstat)

 #------------------------------------------------------------------------------------
 
 ### ------------------------------------------- ###
 ###         Demonstrate generate Meta           ###
 ### ------------------------------------------- ###
 
# multiple of 12 for number of studies -- unbalanced kj
rm(list = ls())
source("SimFunctions/functionsdatagen.R")
dat_kjN <- readRDS("SimFunctions/dat_kjN_mathdat.rds")
set.seed(6535566)
dat_kjN_samp <- n_ES_empirical(dat_kjN, J = 24, with_replacement = TRUE)
design_matrix_ex <- design_matrix(C = 4,
                                  J = 24,
                                  bal = "balanced_j",
                                  k_j = dat_kjN_samp$kj)
design_matrix_ex

mu_vector <- mu_values(
  J = 24,
  tau_sq = .05^2,
  omega_sq = .05^2,
  rho = .5,
  P = .9,
  k_j = dat_kjN_samp$kj,
  N = dat_kjN_samp$N,
  f_c_val = "P5",
  bal = "balanced_j"
)

meta_dat2 <- generate_meta(
  J = 24,
  tau_sq = .05^2,
  omega_sq = .05^2,
  bal = "balanced_j",
  mu_vector = mu_vector,
  rho = .5,
  sample_sizes = dat_kjN_samp$N,
  k_j = dat_kjN_samp$kj,
  f_c_val = "P5",
  return_study_params = FALSE,
  seed = NULL
)


head(meta_dat2)
meta_dat2 |> 
  group_by(category) |>  
  tally()
meta_dat2 |> 
  select(studyid, category) |>  
  distinct()  |> 
  group_by(category)  |>  
  tally()

 #------------------------------------------
 # unbalanced J -- four category , 24 studies
rm(list = ls())
source("SimFunctions/functionsdatagen.R")
dat_kjN <- readRDS("SimFunctions/dat_kjN_mathdat.rds")
set.seed(65436)
dat_kjN_samp <- n_ES_empirical(dat_kjN, J = 24, with_replacement = TRUE)

mu_vector <- mu_values(
  J = 24,
  tau_sq = .05^2,
  omega_sq = .05^2,
  rho = .5,
  P = .9,
  k_j = dat_kjN_samp$kj,
  N = dat_kjN_samp$N,
  f_c_val = "P5",
  bal = "unbalanced_j"
)

meta_dat3 <- generate_meta(
  J = 24,
  tau_sq = .05^2,
  omega_sq = .05^2,
  bal = "unbalanced_j",
  mu_vector = mu_vector,
  rho = .5,
  sample_sizes = dat_kjN_samp$N,
  k_j = dat_kjN_samp$kj,
  f_c_val = "P5",
  return_study_params = FALSE,
  seed = NULL
)
meta_dat3 |> 
  group_by(category) |>  
  tally()

meta_dat3 |> 
  select(studyid, category) |>  
  distinct()  |> 
  group_by(category)  |>  
  tally()

#------------------------------------------
# unbalanced J -- three category
rm(list = ls())
source("SimFunctions/functionsdatagen.R")
dat_kjN <- readRDS("SimFunctions/dat_kjN_mathdat.rds")
dat_kjN_samp <- n_ES_empirical(dat_kjN, J = 24, with_replacement = TRUE)
mu_vector <- mu_values(
  J = 24,
  tau_sq = 0.40^2,
  omega_sq = 0.20^2,
  rho = .5,
  P = .9,
  k_j = dat_kjN_samp$kj,
  N = dat_kjN_samp$N,
  f_c_val = "P4",
  #sigma_j_sq = NA,
  bal = "unbalanced_j"
)

meta_dat3 <- generate_meta(
  J = 24,
  tau_sq = 0.40^2,
  omega_sq = 0.20^2,
  bal = "unbalanced_j",
  mu_vector = mu_vector,
  rho = .5,
  sample_sizes = dat_kjN_samp$N,
  k_j = dat_kjN_samp$kj,
  f_c_val = "P4",
  return_study_params = FALSE,
  seed = NULL
)
meta_dat3  |>  
  group_by(category)  |>  
  tally()

meta_dat3  |>  
  select(studyid, category) |> 
  distinct()  |> 
  group_by(category) |>
  tally()

# just checking if I recapture the tau^2 and omega^2 below
meta_dat3 <- meta_dat3 |>
  group_by(studyid) |>
  mutate(var_g_j = mean(var_g, na.rm = TRUE)) |>
  ungroup()


V_list2 <-
  vcalc(
    vi = meta_dat3$var_g_j,
    cluster = meta_dat3$studyid,
    rho = 0.5,
    obs = meta_dat3$esid
    
    
  )

res_comp <- rma.mv(
  g ~ 0 + category,
  V = V_list2,
  random = ~ 1 | studyid / esid,
  data = meta_dat3,
  test = "t",
  sparse = TRUE,
  verbose = FALSE
)
res_comp
#------------------------------------------

# unbalanced J -- two category
set.seed(21220251)
rm(list = ls())
source("SimFunctions/functionsdatagen.R")
dat_kjN <- readRDS("SimFunctions/dat_kjN_mathdat.rds")
dat_kjN_samp <- n_ES_empirical(dat_kjN, J = 12, with_replacement = TRUE)


design_matrix_ex2 <- design_matrix(C = 2,
                                   J = 12,
                                   bal = "unbalanced_j",
                                   k_j = dat_kjN_samp$kj)
design_matrix_ex2

mu_vector <- mu_values(
  J = 12,
  tau_sq = 05^2,
  omega_sq = .05^2,
  rho = .5,
  P = .9,
  k_j = dat_kjN_samp$kj,
  N = dat_kjN_samp$N,
  f_c_val = "P1",
  #sigma_j_sq = NA,
  bal = "unbalanced_j"
)

meta_dat4 <- generate_meta(
  J = 12,
  tau_sq = .05^2,
  omega_sq = .05^2,
  bal = "unbalanced_j",
  mu_vector = mu_vector,
  rho = .5,
  sample_sizes = dat_kjN_samp$N,
  k_j = dat_kjN_samp$kj,
  f_c_val = "P1",
  return_study_params = FALSE,
  seed = NULL
)

head(meta_dat4)
meta_dat4 |> 
  group_by(category)  |>  
  tally()

meta_dat4  |>  
  select(studyid, category) |> 
  distinct()  |> 
  group_by(category) |> 
  tally()
#------------------------------------------
# sigma_j_q instead of N
set.seed(21220252)
rm(list = ls())
source("SimFunctions/functionsdatagen.R")
dat_kjN <- readRDS("SimFunctions/dat_kjN_mathdat.rds")
dat_kjN_samp <- n_ES_empirical(dat_kjN, J = 12, with_replacement = TRUE)

design_matrix_ex2 <- design_matrix(C = 2,
                                   J = 12,
                                   bal = "unbalanced_j",
                                   k_j = dat_kjN_samp$kj)
design_matrix_ex2

mu_vector <- mu_values(
  J = 12,
  tau_sq = 05^2,
  omega_sq = .05^2,
  rho = .5,
  P = .9,
  k_j = dat_kjN_samp$kj,
  N = dat_kjN_samp$N,
  f_c_val = "P1",
  sigma_j_sq = dat_kjN_samp$se_avg^2,
  bal = "unbalanced_j"
)

meta_dat5 <- generate_meta(
  J = 12,
  tau_sq = .05^2,
  omega_sq = .05^2,
  bal = "unbalanced_j",
  mu_vector = mu_vector,
  rho = .5,
  sample_sizes = dat_kjN_samp$N,
  k_j = dat_kjN_samp$kj,
  f_c_val = "P1",
  return_study_params = FALSE,
  seed = NULL
)

head(meta_dat5)
meta_dat5 |> 
  group_by(category) |> 
  tally()

meta_dat5 |>  
  select(studyid, category) |> 
  distinct()  |> 
  group_by(category) |> 
  tally()

 #------------------------------------------------------------------------------------ 
 ### ------------------------------------------- ###
 ###         Test  Estimate                      ###
 ### ------------------------------------------- ### 
rm(list = ls())
source("SimFunctions/functionsdatagen.R")


set.seed(1232024)

dat_kjN <- readRDS("SimFunctions/dat_kjN_mathdat.rds")
dat_kjN_samp <- n_ES_empirical(dat_kjN, J = 12, with_replacement = TRUE)

mu_vector <- mu_values(
  J = 12,
  tau_sq = 05^2,
  omega_sq = .05^2,
  rho = .5,
  P = .9,
  k_j = dat_kjN_samp$kj,
  N = dat_kjN_samp$N,
  f_c_val = "P5",
  sigma_j_sq = dat_kjN_samp$se_avg^2,
  bal = "balanced_j"
)




meta_dat <- generate_meta(
  J = 12,
  tau_sq = .05^2,
  omega_sq = .05^2,
  bal = "balanced_j",
  mu_vector = mu_vector,
  rho = .5,
  sample_sizes = dat_kjN_samp$N,
  k_j = dat_kjN_samp$kj,
  f_c_val = "P5",
  return_study_params = FALSE,
  seed = NULL
)


test_output <- estimate_model(
  #formula= g ~ 0 + category,
  moderator_val = meta_dat$category,
  cluster_id = meta_dat$studyid,
  delta = meta_dat$g,
  delta_var = meta_dat$var_g,
  es_id = meta_dat$esid,
  r = 0.5,
  smooth_vi = TRUE,
  control_list = list()
)

 
 #smooth variances
 meta_dat <- meta_dat |> 
   group_by(studyid) |>
   mutate(var_g_j = mean(var_g, na.rm = TRUE)) |>
   ungroup()
 
 
 V_list2 <- 
   vcalc(
     vi = meta_dat$var_g_j,
     cluster = meta_dat$studyid,
     rho = 0.5,
     obs = meta_dat$esid

     
   )
 
 res_comp <- rma.mv(g ~ 0 + category,
   V = V_list2, 
   random = ~ 1 | studyid / esid,
   data = meta_dat,
   test = "t",
   sparse = TRUE,
   verbose = FALSE
 )
 
 coef_RVE2 <-  metafor::robust(
   res_comp, # estimation model above
   cluster = studyid, # define clusters
   clubSandwich = TRUE # use CR2 adjustment
 )
 
 wald_test_results2 <- Wald_test((res_comp), 
                                constraints = constrain_equal(1:4), 
                                vcov =  "CR2")
 

 
 rownames(coef_RVE2$beta) <- NULL
 metafor_output <- coef_RVE2$beta[1:4, 1]

 myfunc_output <- unlist(test_output$est) 

 
 all.equal(myfunc_output,metafor_output)
 all.equal(test_output$est_var, list(coef_RVE2$se^2))
 
 
 #----------------------------------------------------
 # test if it uses rma.uni() when only kj = 1
 rm(list = ls())
 source("SimFunctions/functionsdatagen.R")
 
 
 set.seed(1232024)
 
 dat_kjN <- readRDS("SimFunctions/dat_kjN_mathdat.rds")
 dat_kjN_samp <- n_ES_empirical(dat_kjN, J = 12, with_replacement = TRUE)
 
 mu_vector <- mu_values(
   J = 12,
   tau_sq = 05^2,
   omega_sq = .05^2,
   rho = .5,
   P = .9,
   k_j = dat_kjN_samp$kj,
   N = dat_kjN_samp$N,
   f_c_val = "P5",
   sigma_j_sq = dat_kjN_samp$se_avg^2,
   bal = "balanced_j"
 )
 
 meta_dat <- generate_meta(
   J = 12,
   tau_sq = .05^2,
   omega_sq = .05^2,
   bal = "balanced_j",
   mu_vector = mu_vector,
   rho = .5,
   sample_sizes = dat_kjN_samp$N,
   k_j = dat_kjN_samp$kj,
   f_c_val = "P5",
   return_study_params = FALSE,
   seed = NULL
 )
 
 meta_dat_uni <- meta_dat |>
   group_by(studyid) |>
   slice(1:1) |>
   ungroup()
 
 
 test_output2 <- estimate_model(
   moderator_val = meta_dat_uni$category,
   cluster_id = meta_dat_uni$studyid,
   delta = meta_dat_uni$g,
   delta_var = meta_dat_uni$var_g,
   es_id = meta_dat_uni$esid,
   r = 0.5,
   smooth_vi = TRUE,
   control_list = list()
 )
 
 #----------------------------------------------------
 # test estimate that will result in a non-convergence
 
 
 rm(list=ls())
 source("SimFunctions/functionsdatagen.R")
 set.seed(522719220)
 
 
 meta_dat <- data.frame(
   studyid = 1:10,
   g  = rep(.2, 10),
   var_g = rep(.02, 10),
   category = c(rep("A", 5),rep("B", 5)),
   esid = 1:10
 ) 
 
 
 test_output <- estimate_model( 
   data =  meta_dat,
   moderator_val = meta_dat$category,
   cluster_id = meta_dat$studyid,
   delta = meta_dat$g, 
   delta_var = meta_dat$var_g,es_id = meta_dat$esid, 
   r= 0.5, 
   smooth_vi = TRUE, 
   control_list = list()
 )
 

 all.equal(test_output$dat[[1]], meta_dat)

 
 #------------------------------------------------------------------------------------ 
 ### ------------------------------------------- ###
 ###         Demonstrate  run_sim                ###
 ### ------------------------------------------- ### 
 rm(list=ls())
 source("SimFunctions/functionsdatagen.R")



 dat_kjN <- readRDS("SimFunctions/dat_kjN_mathdat.rds")
 
 
 tm1 <- system.time(test <- run_sim(iterations = 1, 
                                    J = 60, tau_sq = (0.40)^2, 
                                    omega_sq = (0.20)^2, 
                                    bal = "balanced_j", 
                                    rho = .8, P = .9, 
                                    f_c_val = "P8",
                                    sigma_j_sq_inc = FALSE,
                                    pilot_data = dat_kjN, 
                                    return_study_params = FALSE,
                                    seed = NULL,
                                    summarize_results = FALSE))
 
 tm2 <- system.time(test2 <- run_sim(iterations = 3, 
                                      J = 60, tau_sq = (0.40)^2, 
                                      omega_sq = (0.20)^2, 
                                      bal = "balanced_j", 
                                      rho = .8, P = .9, 
                                      f_c_val = "P8",
                                      sigma_j_sq_inc = FALSE,
                                      pilot_data = dat_kjN, 
                                      return_study_params = FALSE,
                                      seed = NULL,
                                      summarize_results = FALSE))
 

 
 ### Notes
 #number of conditions
 5*2*2*2*6*8*2 
 #3840
 
 #time in hours for run (with J=72) (one core/4.80 GHz) * number reps (2500) * conditions  node hours (48 cores per node)
 ((tm1[[3]]*2500*3840)/48)/60^2
 
 

 