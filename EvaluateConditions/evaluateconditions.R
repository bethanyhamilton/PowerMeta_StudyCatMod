source("SimFunctions/functionsdatagen.R")



design_factors_bal <- list(
  C = c(2, 3, 4),
  J = c(12, 24, 36, 48, 60, 72),
  tau_sq = c(0.05, 0.20, 0.40)^2, 
  omega_sq = c(0.05, 0.10, 0.20)^2,
  rho = c(.2, .5, .8),
  P = c(0.05, 0.2, 0.4, 0.6, 0.8, 0.9)
)


expand.grid(design_factors_bal)





run_power <- function(C, 
                      J, 
                      tau_sq, 
                      omega_sq, 
                      rho, 
                      k_j = 3, ####CHANGE THIS LATER
                      n_j = 30, ####CHANGE THIS LATER
                      mu_values
                      ){
  
 dat_app <-  dat_approx(C = C, J = J, tau_sq = tau_sq, 
                        omega_sq = omega_sq, rho = rho, 
                        k_j = k_j, n_j = n_j )
 
 
 power <-  power_CHE_RVE_study_cat(data = dat_app, 
                         moderator_val  = cat, 
                         cluster_id = studyid, 
                         sigma_j_sq_val = sigma_j_sq,
                         rho_val = rho,
                         omega_sq_val = omega_sq,
                         tau_sq_val = tau_sq,
                         mu = mu_values, 
                         alpha = .05)
 


 return(power)
  
}

mu_values <- function(C, 
                      J, 
                      tau_sq, 
                      omega_sq, 
                      rho, 
                      P, 
                      k_j = 3, ####CHANGE THIS LATER
                      n_j = 30, ####CHANGE THIS LATER
                      f_c_val ) {
  
  dat_app <-  dat_approx(C = C, J = J, tau_sq = tau_sq, 
                         omega_sq = omega_sq, rho = rho, 
                         k_j = k_j, n_j = n_j )
  
  dfs <- multiple_categories(data = dat_app, 
                             moderator_val  = cat, 
                             cluster_id = studyid, 
                             sigma_j_sq_val = sigma_j_sq,
                             rho_val = rho,
                             omega_sq_val = omega_sq,
                             tau_sq_val = tau_sq)
  
  df_num <- dfs$df_num[1]
  df_den <- dfs$nu_D[1] - dfs$df_num[1] + 1
  
  lambda <- find_lambda(d1= df_num, 
                        d2 = df_den, 
                        x = qf(1-.05, df_num,df_den ), 
                        area = 1 - P, 
                        interval = c(0, 100), 
                        tol = 0.0001)
    
  zeta_val <- zeta(pattern = f_c_val, lambda = lambda, weights =dfs$W)
  
  mu_values <- build_mu(pattern = f_c_val, zeta = zeta_val)
  
  return(mu_values)
  
  
}

mu_vector <- mu_values(C = 4, J = 12, tau_sq = .05^2, omega_sq = .05^2, 
                       rho = .5, P = .3, k_j = 3, n_j = 30, f_c_val = f_c(pattern = 4))


run_power(C = 4,
          J = 12,
          tau_sq = .05^2,
          omega_sq = .05^2,
          rho = .5,
          k_j = 3,
          ####CHANGE THIS LATER
          n_j = 30,
          ####CHANGE THIS LATER
          mu_values= mu_vector)

#apply function mu to row of the combinations of conditions. what sort of mu. 

# ## 3 categories, J_c = 3
# data_ex <- tibble(studyid = c(1:9),
#                   k_j = rep(3, 9),
#                   n_j = rep(30, 9),
#                   sigma_j = sqrt(4 / n_j), 
#                   omega = rep(.10, 9),
#                   rho = rep(.5, 9),
#                   tau = rep(.20, 9),
#                   three_cat =c(rep("a",3), rep("b",3), rep("c",3))
# )
# 
# 
# df2 <- multiple_categories(dat = data_ex, moderator = data_ex$three_cat, cluster = data_ex$studyid, sigma_j = data_ex$sigma_j , rho = data_ex$rho, omega_sq = data_ex$omega, tau_sq = data_ex$tau)
# 
# df_num <- df2$df_num[1]
# df_den <- df2$nu_D[1] - df2$df_num[1] + 1
# 
# f_crit <- qf(1-.05, df_num,df_den )
# 
# lambda <- find_lambda(d1= df_num, d2 = df_den, x = f_crit, area = 1 - .9, interval = c(0, 100), tol = 0.0001)
# 
# f_c <- c(0, 1, 2)
# #f_c <- c(0, 0, 1)
# #f_c <- c(0, 1, 1)
# 
# zeta_val <- zeta(pattern = f_c, lambda = lambda, weights =df2$W )
# 
# mu_values <- build_mu(pattern = f_c, zeta = zeta_val)
# 
# # test to get the same. 
# power_CHE_RVE_study_cat(dat = data_ex, moderator = data_ex$three_cat, cluster = data_ex$studyid, sigma_j = data_ex$sigma_j , rho = data_ex$rho, omega_sq = data_ex$omega, tau_sq = data_ex$tau, mu = mu_values, alpha = .05)
# 
# # 


# 
## test find zeta
#data_ex <- tibble(studyid = c(1:40),
#                  k_j = rep(10, 40),
#                  n_j = rep(100, 40),
#                  sigma_j = sqrt(4 / n_j), 
#                  omega = rep(.10, 40),
#                  rho = rep(.5, 40),
#                  tau = rep(.20, 40),
#                  #   category = sample(letters[1:C], size = 40, replace = TRUE, prob = c(.25, .25,.25, .25))
#                  two_cat =c(rep("a",20), rep("b",20))
#)

#data_ex$four_cat <- c(rep("a",10), rep("b",10), rep("c",10), rep("d",10))


#df2 <- multiple_categories(dat = data_ex, moderator = data_ex$four_cat, cluster = data_ex$studyid, sigma_j = data_ex$sigma_j, rho = data_ex$rho,omega_sq = data_ex$omega, tau_sq = data_ex$tau)


#df_num <- df2$df_num[1]
#df_den <- df2$nu_D[1] - df2$df_num[1] + 1

#weights <- df2$W

#fp <- c(0, 1, 2, 3)

#zeta_val <- zeta(pattern = fp, lambda = 1.562499, weights =weights )




#build_mu(pattern = fp, zeta = zeta)

