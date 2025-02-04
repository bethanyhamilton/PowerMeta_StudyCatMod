#rm(list=ls())


library(tidyverse)
library(metafor)
library(clubSandwich)
library(mvtnorm)

#source("SimFunctions/functionsdatagen.R")

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


# create mu vector
# dat_approx(C = 3, J = 12, tau_sq = .05^2, 
#            omega_sq = .05^2, rho = .5, 
#            k_j = 3,
#            N = 30, 
#         #   sigma_j_sq = NA,
#            bal = "balanced" 
# )

rm(list=ls())

source("SimFunctions/functionsdatagen.R")

## make my vector
mu_vector <- mu_values(J = 12, tau_sq = .05^2, 
                       omega_sq = .05^2,
                       rho = .5, P = .9, 
                       k_j = 3, N = 30, 
                       f_c_val = 5, 
                       #sigma_j_sq = NA,
                       bal ="balanced" )
mu_vector


## see if I get same power as I used to generate mu when I plug in mu values
run_power(C = 4,
          J = 12,
          tau_sq = .05^2,
          omega_sq = .05^2,
          rho = .5,
          k_j = 3,
          N = 30,
          bal = "balanced",
          mu_values= mu_vector)




#------------------------------------------------------------------------------------

### ------------------------------------------- ###
###         Test generate Meta and Df           ###
### ------------------------------------------- ###
rm(list=ls())

source("SimFunctions/functionsdatagen.R")


sample_empirical_dat <- tibble(N = rep(200, 12 ), k_j  = c(3,4,3,5,3,3,3,3,3,3,3,3) )
meta_dat <- generate_meta(J = 12, tau_sq = .05^2, 
                          omega_sq = .05^2, bal = "balanced", C = 4,
                                      # cor_mu, cor_sd,
                                      rho = .5, P = .9, sample_sizes = sample_empirical_dat$N, k_j = sample_empirical_dat$k_j,
                                      ### added k_j here for now, but it will probably
                                      # be part of sample_sizes so should remove later.
                                       f_c_val = 5,
                                      return_study_params = FALSE,
                                      seed = NULL)





### clubSandwich df
V <- with(meta_dat, 
          clubSandwich::impute_covariance_matrix(vi = var_g, cluster = studyid, r = .5, smooth_vi = TRUE))

fit <- rma.mv(g ~ -1 + category, 
                      V = V, random = ~ 1 | studyid/esid, 
                      data = meta_dat, test = "t", method = "REML")

V_sep <- vcovCR(fit, cluster = meta_dat$studyid, type = "CR2")

clubsandwichdf<- Wald_test((fit), constraints = constrain_equal(1:4), vcov = V_sep)
#clubsandwichdf$df_denom
 

### my function df 
df  <- multiple_categories(data = meta_dat,
                    moderator_val = category,
                    cluster_id = studyid,
                    sigma_j_sq_val = var_g,
                    rho_val = .5,
                    omega_sq_val = fit$sigma2[2],
                    tau_sq_val = fit$sigma2[1]

                      )

df_Z <- as.numeric(df[1,5])
q <- as.numeric(df[1,6])
C_sep <- constrain_equal(1:4, coefs = coef(fit))
mu_hat <- coef(fit)
VR <- vcovCR(fit, type = "CR2") %>% as.matrix()


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
 
 
 
 
 
 #------------------------------------------------------------------------------------ 
 ### ------------------------------------------- ###
 ###         Test  Estimate                      ###
 ### ------------------------------------------- ### 
 rm(list=ls())
 source("SimFunctions/functionsdatagen.R")
 

 set.seed(1232024)
 test_dat <- tibble(N = rep(200, 12 ), k_j  = c(3,4,3,5,3,3,3,3,3,3,3,3) )
 meta_dat <- generate_meta(J = 12, tau_sq = .05^2, omega_sq = .05^2, bal = "balanced", C = 4,
                          rho = .5, P = .9, sample_sizes = test_dat$N, k_j = test_dat$k_j,
                          f_c_val = 5, return_study_params = FALSE,
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
 
 
 
 all.equal(test_output$est, list(coef_RVE2$beta))
 all.equal(test_output$est_var, list(coef_RVE2$se^2))
 

 
 
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
# data_ex <- tibble(studyid = c(1:40),
#                  k_j = rep(10, 40),
#                  n_j = rep(100, 40),
#                  sigma_j = sqrt(4 / n_j),
#                  omega = rep(.10, 40),
#                  rho = rep(.5, 40),
#                  tau = rep(.20, 40),
#                  #   category = sample(letters[1:C], size = 40, replace = TRUE, prob = c(.25, .25,.25, .25))
#                  two_cat =c(rep("a",20), rep("b",20))
# )
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

