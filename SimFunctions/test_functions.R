library(tidyverse)
library(metafor)
library(clubSandwich)
library(mvtnorm)


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


# run_power(C = 4,
#           J = 12,
#           tau_sq = .05^2,
#           omega_sq = .05^2,
#           rho = .5,
#           k_j = 3,
#           N = 30,
#           bal = "balanced_j",
#           mu= mu_vector)

## need to construct the mu to find power, so this function now gets the Power and 
# f_c_val (pattern) of mu values to make the mu vector. Then it uses the constructed
# mu vector to find power. It matches up. I needed to do this because later I have a vector
# of k_j and N values from different sampling methods and I needed the mu vector to change with those
# values. 


run_power(C = 4,
          J = 12,
          tau_sq = .05^2,
          omega_sq = .05^2,
          rho = .5,
          k_j = 3,
          N = 30,
          bal = "balanced_j",
          P = .9, 
          f_c_val = "P5")



### --------------------------------------------------------- ###
### Demonstrate Power Approximation given a set of conditions ###
### --------------------------------------------------------- ###
rm(list=ls())
source("SimFunctions/functionsdatagen.R")
set.seed(2202025)

dat_kjN <- readRDS("SimFunctions/dat_kjN_erikadat.rds")
dat_kjN_samp <- n_ES_empirical(dat_kjN, J = 24, with_replacement = TRUE)
shape_rate <- MASS::fitdistr(dat_kjN$N, "gamma")

# test <- power_approximation(C = 4, 
#                     J = 24, 
#                     tau_sq = .40^2, 
#                     omega_sq= .10^2, 
#                     rho = 0.8, 
#                     sigma_j_sq = NULL,
#                     N_mean = mean(dat_kjN$N), 
#                     k_mean = mean(dat_kjN$kj),
#                     N_dist = shape_rate, 
#                     pilot_data = dat_kjN_samp, 
#                     iterations = 10,
#                     sample_size_method = c("balanced","stylized","empirical"),
#                     mu_vec = mu_vector,
#                     bal = "balanced_j",
#                     seed = NULL)
# test


test <- power_approximation(C = 4, 
                            J = 24, 
                            tau_sq = .40^2, 
                            omega_sq= .10^2, 
                            rho = 0.8, 
                            sigma_j_sq = NULL,
                            N_mean = mean(dat_kjN$N), 
                            k_mean = mean(dat_kjN$kj),
                            N_dist = shape_rate, 
                            pilot_data = dat_kjN_samp, 
                            iterations = 10,
                            sample_size_method = c("balanced","stylized","empirical"),
                            #mu_vec = mu_vector,
                            P = .5,
                           # f_c_val = 5,
                           f_c_val = "P5",
                            bal = "balanced_j",
                            seed = NULL)
test



#------------------------------------------------------------------------------------

### ------------------------------------------- ###
###         Test generate Meta and Df           ###
### ------------------------------------------- ###
rm(list=ls())
source("SimFunctions/functionsdatagen.R")


sample_empirical_dat <- tibble(N = rep(200, 12 ), k_j  = c(3,4,3,5,3,3,3,3,3,3,3,3) )
meta_dat <- generate_meta(J = 12, tau_sq = .05^2, 
                          omega_sq = .05^2, 
                          bal = "balanced_j", 
                          #C = 4,
                                      rho = .5, P = .9, sample_sizes = sample_empirical_dat$N, 
                                      k_j = sample_empirical_dat$k_j,
                                      sigma_j_sq = NULL,
                                      f_c_val = "P5",
                                      return_study_params = FALSE,
                                      seed = NULL)




### clubSandwich df
V <- with(meta_dat, 
          clubSandwich::impute_covariance_matrix(vi = var_g, cluster = studyid, r = .5, smooth_vi = TRUE))

fit <- rma.mv(g ~ -1 + category, 
                      V = V, random = ~ 1 | studyid/esid, sparse = TRUE,
                      data = meta_dat, test = "t", method = "REML")


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

df_Z <- df_val(E_Vr = multiplecat$E_Vr, nu_Vr = multiplecat$nu_Vr, W = multiplecat$W)[[2]]
q <- df_val(E_Vr = multiplecat$E_Vr, nu_Vr = multiplecat$nu_Vr, W = multiplecat$W)[[1]]

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
 rm(list=ls())
 source("SimFunctions/functionsdatagen.R")
 dat_kjN <- readRDS("SimFunctions/dat_kjN_erikadat.rds")
 set.seed(6535566)
 


 dat_kjN_samp <- n_ES_empirical(dat_kjN, J = 24, with_replacement = TRUE)
 
 design_matrix_ex <- design_matrix(C = 4, J = 24, bal = "balanced_j",  k_j = dat_kjN_samp$kj )
 design_matrix_ex
 
 meta_dat2 <- generate_meta(J = 24, tau_sq = .05^2, 
                            omega_sq = .05^2, 
                            bal = "balanced_j", 
                            #C = 4,
                            rho = .5, P = .9, sample_sizes = dat_kjN_samp$N, 
                            k_j = dat_kjN_samp$kj,
                            sigma_j_sq = NULL,
                            f_c_val = "P5",
                            return_study_params = FALSE,
                            seed = NULL)
 
 
 head(meta_dat2)
 meta_dat2 |> group_by(category) |>  tally()
 meta_dat2|> select(studyid, category) |>  distinct()  |> group_by(category)  |>  tally()
 
 #------------------------------------------
 # unbalanced J -- four category , 24 studies
 rm(list=ls())
 source("SimFunctions/functionsdatagen.R")
 dat_kjN <- readRDS("SimFunctions/dat_kjN_erikadat.rds")
 set.seed(65436)
 dat_kjN_samp <- n_ES_empirical(dat_kjN, J = 24, with_replacement = TRUE)
 meta_dat3 <- generate_meta(J = 24, tau_sq = .05^2, 
                            omega_sq = .05^2, 
                            bal = "unbalanced_j", 
                            #C = 4,
                     
                            rho = .5, P = .9, sample_sizes = dat_kjN_samp$N, 
                            k_j = dat_kjN_samp$kj,
                            sigma_j_sq = NULL,
                            f_c_val = "P5",
                            return_study_params = FALSE,
                            seed = NULL)
 meta_dat3 |> group_by(category) |>  tally()
 
 meta_dat3|> select(studyid, category) |>  distinct()  |> group_by(category)  |>  tally()
 
 #------------------------------------------
 # unbalanced J -- three category
 
 mu_vector <- mu_values(J = 24, tau_sq = 0.40^2, 
                        omega_sq = 0.20^2,
                        rho = .5, P = .9, 
                        k_j = dat_kjN_samp$kj, N = dat_kjN_samp$N, 
                        f_c_val = "P4", 
                        #sigma_j_sq = NA,
                        bal ="unbalanced_j" )
 mu_vector
 
 meta_dat3 <- generate_meta(J = 24, tau_sq = 0.40^2, 
                            omega_sq = 0.20^2, 
                            bal = "unbalanced_j", 
                            #C = 3,
                            rho = .5, P = .9, sample_sizes = dat_kjN_samp$N, 
                            k_j = dat_kjN_samp$kj,
                            sigma_j_sq = NULL,
                            f_c_val = "P4",
                            return_study_params = FALSE,
                            seed = NULL)
 meta_dat3  |>  group_by(category)  |>  tally()
 
 meta_dat3  |>  select(studyid, category) |>  distinct()  |> group_by(category) |> tally()
 
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
 
 res_comp <- rma.mv(g ~ 0 + category,
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
 
 dat_kjN_samp2 <- n_ES_empirical(dat_kjN_samp, J = 12, with_replacement = TRUE)

 design_matrix_ex2 <- design_matrix(C = 2, J = 12, bal = "unbalanced_j",  k_j = dat_kjN_samp2$kj )
 design_matrix_ex2
 
 
 meta_dat4 <- generate_meta(J = 12, tau_sq = .05^2, 
                            omega_sq = .05^2, 
                            bal = "unbalanced_j",
                            #C = 2,
                            rho = .5, P = .9, sample_sizes = dat_kjN_samp2$N, 
                            k_j = dat_kjN_samp2$kj,
                            sigma_j_sq = NULL,
                            f_c_val = "P1",
                            return_study_params = FALSE,
                            seed = NULL)
 
 head(meta_dat4)
 meta_dat4 |> group_by(category)  |>  tally()
 
 meta_dat4  |>  select(studyid, category)|> distinct()  |> group_by(category)|> tally()
 
 
 #------------------------------------------
 # sigma_j_q instead of N
 set.seed(21220252)
 dat_kjN_samp2 <- n_ES_empirical(dat_kjN_samp, J = 12, with_replacement = TRUE)
 
 design_matrix_ex2 <- design_matrix(C = 2, J = 12, bal = "unbalanced_j",  k_j = dat_kjN_samp2$kj )
 design_matrix_ex2
 
 
 meta_dat5 <- generate_meta(J = 12, tau_sq = .05^2, 
                            omega_sq = .05^2, 
                            bal = "unbalanced_j", 
                            #C = 2,
                            rho = .5, P = .9, sample_sizes = dat_kjN_samp2$N, 
                            k_j = dat_kjN_samp2$kj,
                            sigma_j_sq = dat_kjN_samp2$se_avg^2,
                            f_c_val = "P1",
                            return_study_params = FALSE,
                            seed = NULL)
 
 head(meta_dat5)
 meta_dat5 |> group_by(category)  |>  tally()
 
 meta_dat5  |>  select(studyid, category)|> distinct()  |> group_by(category)|> tally()
 #------------------------------------------
 
 #### this seed resulted in an error. Some combination of k_j and N result in 
 #### an error in cov_mat object
 
 ### seed with error...the study with 14 effects
 ## k_j=14, N= 6, Sigma = .5
 #Sigma <- Sigma + diag(1 - Sigma, nrow = k_j)
 #cov_mat <- as.matrix(rWishart(n = 1, df = N - 2, Sigma = Sigma)[,,1])
 
 #Solution: df >= dimension
 #NEED TO PREPROCESS DATA SO  N >= kj + 2

 rm(list=ls())
 
 ## df greater kj 
 source("SimFunctions/functionsdatagen.R")
 dat_kjN <- readRDS("SimFunctions/dat_kjN_erikadat.rds")
 set.seed(2122025)
 dat_kjN_samp <- n_ES_empirical(dat_kjN, J = 24, with_replacement = TRUE)
 
 dat_kjN_samp <- dat_kjN_samp |>  filter(kj == 12)
 meta_dat2 <- generate_meta(J = 24, tau_sq = .05^2, 
                            omega_sq = .05^2, 
                            bal = "balanced_j", 
                            #C = 4,
                            rho = .5, P = .9, sample_sizes = dat_kjN_samp$N, 
                            k_j = dat_kjN_samp$kj,
                            sigma_j_sq = NULL,
                            f_c_val = "P5",
                            return_study_params = FALSE,
                            seed = NULL)
 

 head(meta_dat2)
 

 #------------------------------------------------------------------------------------ 
 ### ------------------------------------------- ###
 ###         Test  Estimate                      ###
 ### ------------------------------------------- ### 
 rm(list=ls())
 source("SimFunctions/functionsdatagen.R")
 

 set.seed(1232024)
 
 test_dat <- tibble(N = rep(200, 12 ), k_j  = c(3,4,3,5,3,3,3,3,3,3,3,3) )
 meta_dat <- generate_meta(J = 12, tau_sq = .05^2, omega_sq = .05^2, bal = "balanced_j", 
                           #C = 4,
                          rho = .5, P = .9, sample_sizes = test_dat$N, k_j = test_dat$k_j,
                          f_c_val = "P5", return_study_params = FALSE,
                                       seed = NULL)

 # test_output <- estimate_model(data= meta_dat, formula= g ~ 0 +category, 
 #                               moderator_val = category,cluster_id = studyid,
 #                               delta = g, delta_var = var_g,es_id = esid, r= 0.5, 
 #                       smooth_vi = TRUE, control_list = list()
 # )
 
 test_output <- estimate_model( 
   #formula= g ~ 0 + category, 
                               moderator_val = meta_dat$category,
                               cluster_id = meta_dat$studyid,
                               delta = meta_dat$g, 
                               delta_var = meta_dat$var_g,es_id = meta_dat$esid, 
                               r= 0.5, 
                               smooth_vi = TRUE, control_list = list()
 )

 
# meta_dat$studyid <- as.factor(meta_dat$studyid)
 
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
 
 coef_RVE2 <-  robust(
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
 
 
 #------------------------------------------------------------------------------------ 
 ### ------------------------------------------- ###
 ###         Demonstrate  run_sim                ###
 ### ------------------------------------------- ### 
 rm(list=ls())
 source("SimFunctions/functionsdatagen.R")
 #dat_kjN <- readRDS("SimFunctions/dat_kjN.rds")
 
 
 ## need to test with different empirical data set that has SE
 ## also should capture the mu_values from the data generation. 

 dat_kjN <- readRDS("SimFunctions/dat_kjN_Diet_dat.rds")
 
 test <- run_sim(iterations = 1, 
                 J = 24, tau_sq = .05^2, 
                 omega_sq = .05^2, 
                 bal = "balanced_j", 
                 rho = .5, P = .9, 
                 f_c_val = "P5",
                 sigma_j_sq_inc = FALSE,
                 pilot_data = dat_kjN, 
                 return_study_params = FALSE,
                 seed = NULL,
                 summarize_results = FALSE)
 
 tm1 <- system.time(test <- run_sim(iterations = 1, 
                                    J = 72, tau_sq = .05^2, 
                                    omega_sq = .05^2, 
                                    bal = "balanced_j", 
                                    rho = .5, P = .9, 
                                    f_c_val = "P5",
                                    sigma_j_sq_inc = FALSE,
                                    pilot_data = dat_kjN, 
                                    return_study_params = FALSE,
                                    seed = NULL,
                                    summarize_results = FALSE))
 
 
 tm2 <- system.time(test2 <- run_sim2(iterations = 1, 
                                    J = 72, tau_sq = .05^2, 
                                    omega_sq = .05^2, 
                                    bal = "balanced_j", 
                                    rho = .5, P = .9, 
                                    f_c_val = "P5",
                                    sigma_j_sq_inc = FALSE,
                                    pilot_data = dat_kjN, 
                                    return_study_params = FALSE,
                                    seed = NULL,
                                    summarize_results = FALSE))
 
 
 #number of conditions
 5*2*2*2*6*8*2 
 #3840
 
 #time for run (with J=72) * number reps (2500) * conditions
 tm1[[3]]*2500*3840
 
 tm2[[3]]*2500*3840
 
# df_Z <- multiple_categories(dat = meta_dat, moderator = meta_dat$category, cluster = meta_dat$studyid, fit = fit, d_var = meta_dat$var_g, rho = .5)
# 


# multiple_categories <- function(dat, moderator, cluster, fit, d_var , rho){
#   
#   c <-  length(coef(fit))
#   q <-  c-1
#   
#   dat <- dat %>% 
#     mutate(moderator = moderator,
#            cluster = cluster, 
#            d_var = d_var)
#   
#   
#   study_level_data <-  dat %>% group_by(moderator, cluster) %>% summarise(mean_var = mean(d_var), k_j = n(), .groups = 'drop') 
#   
#   study_level_data$tau_sq <- fit$sigma2[1]
#   study_level_data$omega_sq <- fit$sigma2[2]
#   study_level_data$rho <- rho
#   
#   by_category <-  study_level_data %>% 
#     group_by(moderator) %>% 
#     mutate(
#       w_j = study_level_weights(k_j =k_j ,sigma_j = mean_var, tau_sq = tau_sq, rho = rho, omega_sq = omega_sq )) %>% 
#     summarise(W = sum(w_j),
#               E_Vr = 1/W,
#               nu_Vr = (sum(w_j^2/(W-w_j)^2) - ((2/W)*sum(w_j^3/(W-w_j)^2)) + ((1/W^2)*(sum(w_j^2/(W-w_j)))^2)) ^(-1), .groups = 'drop' )
#   
#   
#   
#   nu_D <- (q*(q+1)) / (2*(sum( (1/by_category$nu_Vr)* (1-(1/(by_category$E_Vr*sum(by_category$W))))^2 )))
#   
#   return(nu_D)
#   
# }
# 
# q = 4-1



 # 
# df_num <- df2$df_num[1]
# df_den <- df2$nu_D[1] - df2$df_num[1] + 1
# # 
#  f_crit <- qf(1-.05, df_num,df_den )
# # 
#  lambda <- find_lambda(d1= df_num, d2 = df_den, x = f_crit, area = 1 - .9, interval = c(0, 100), tol = 0.0001)
#  # 
#  f_c <- c(0, 1, 2)
# # #f_c <- c(0, 0, 1)
# # #f_c <- c(0, 1, 1)
# # 
#  zeta_val <- zeta(pattern = f_c, lambda = lambda, weights =df2$W )
# # 
#  mu_values <- build_mu(pattern = f_c, zeta = zeta_val)
# # 
# # # test to get the same. 
#  power_CHE_RVE_study_cat(dat = data_ex, moderator = data_ex$three_cat, cluster = data_ex$studyid, sigma_j = data_ex$sigma_j , rho = data_ex$rho, omega_sq = data_ex$omega, tau_sq = data_ex$tau, mu = mu_values, alpha = .05)
# # 
# # # 
# 
# 
# # 
# ## test find zeta
#data_ex <- tibble(studyid = c(1:40),
#                 k_j = rep(10, 40),
#                 n_j = rep(100, 40),
#                 sigma_j = sqrt(4 / n_j),
#                 omega = rep(.10, 40),
#                 rho = rep(.5, 40),
#                 tau = rep(.20, 40),
#                 #   category = sample(letters[1:C], size = 40, replace = TRUE, prob = c(.25, .25,.25, .25))
#                 two_cat =c(rep("a",20), rep("b",20))
#)
# 
# data_ex$four_cat <- c(rep("a",10), rep("b",10), rep("c",10), rep("d",10))
# 
# 
# df2 <- multiple_categories(dat = data_ex, moderator = data_ex$four_cat, cluster = data_ex$studyid, sigma_j = data_ex$sigma_j, rho = data_ex$rho,omega_sq = data_ex$omega, tau_sq = data_ex$tau)
# 
# 
# df_num <- df2$df_num[1]
# df_den <- df2$nu_D[1] - df2$df_num[1] + 1

#weights <- df2$W

#fp <- c(0, 1, 2, 3)

#zeta_val <- zeta(pattern = fp, lambda = 1.562499, weights =weights )




#build_mu(pattern = fp, zeta = zeta)

