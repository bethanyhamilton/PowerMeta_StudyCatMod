
library(tidyverse)

# study level weights
study_level_weights <- function(k_j,sigma_j, tau_sq, rho, omega_sq  ){
  weights <- k_j/ (k_j*tau_sq + (k_j-1)*rho*sigma_j + omega_sq + sigma_j)
  
  return(weights)
  
  
}





multiple_categories <- function(dat, moderator, cluster, sigma_j , rho, omega_sq, tau_sq){
  
  
  
  dat <- dat %>% 
    mutate(moderator = moderator,
           cluster = cluster, 
           sigma_j = sigma_j,
           omega_sq = omega_sq,
           tau_sq = tau_sq,
           rho = rho)
  
  c <-  length(unique(dat$moderator))
  q <-  c-1
  
  df_num <- q
  
  study_level_data <-  dat %>% 
    group_by(cluster) %>% 
    mutate(k_j = n()) %>% ungroup()

  
  
  by_category <-  study_level_data %>% 
    group_by(moderator) %>% 
    mutate(
      w_j = study_level_weights(k_j =k_j ,sigma_j = sigma_j, tau_sq = tau_sq, rho = rho, omega_sq = omega_sq )) %>% 
    summarise(W = sum(w_j),
              E_Vr = 1/W,
              nu_Vr = (sum(w_j^2/(W-w_j)^2) - ((2/W)*sum(w_j^3/(W-w_j)^2)) + ((1/W^2)*(sum(w_j^2/(W-w_j)))^2)) ^(-1), .groups = 'drop' )
  
  
  
  by_category$nu_D <- (q*(q+1)) / (2*(sum( (1/by_category$nu_Vr)* (1-(1/(by_category$E_Vr*sum(by_category$W))))^2 )))
  by_category$df_num = df_num
  
  
  return(by_category)
  
}

ncp <- function(weights, mu_p){
  
  mu_wavg = sum(weights*(mu_p))/ sum(weights)
  ncp = sum(weights*(mu_p-mu_wavg)^2)
  return(ncp)
}




power <- function(dat, moderator, cluster, sigma_j , rho, omega_sq, tau_sq, mu, alpha) {
  
  df <- multiple_categories(dat = dat, moderator = moderator, cluster = cluster,sigma_j = sigma_j, rho = rho,omega_sq = omega_sq, tau_sq = tau_sq)
  
  ncp <-  ncp(weights = df$W, mu_p = mu)
  
  df_num <- df$df_num[1]
  df_den <- df$nu_D[1] - df$df_num[1] + 1
  
  
  F_crit <- qf(1-alpha, df_num, df_den)
  power <- 1- pf(F_crit, df_num, df_den, ncp = ncp)
  
  
  power_ncp <- cbind(power, ncp)
  
  return(power_ncp)
  
  
}


# plug into uniroot 
non_central_f_cdf_reverse <- function(lambda, d1, d2, x, area){
  
  pf(x, d1, d2, ncp = lambda) - area
  
  
}

# find lambda

find_lambda <- function(lambda, d1, d2, x, area, interval, tol){
  
 roots <- uniroot(f = non_central_f_cdf_reverse, d1 = d1, d2 = d2, x = x, area = area, interval = interval,  extendInt = "yes", tol = tol)
  
 return(roots$root)
 
}

#find_lambda(d1 = 3, d2 = 23, x = 3.027998, area = 0.8603392, interval = c(0, 10), tol = 0.0001)


zeta <- function(pattern, lambda, weights){
  
  zeta <- sqrt((lambda)/(sum(weights*(pattern - sum(weights*pattern)/sum(weights) )^2)))
  
  return(zeta)
}


build_mu <- function(pattern, zeta){
  
  fp*zeta
  
}


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
