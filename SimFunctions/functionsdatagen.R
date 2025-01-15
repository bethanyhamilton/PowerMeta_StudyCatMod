

# study level weights
study_level_weights <- function(k_j,sigma_j_sq, 
                                tau_sq, rho, omega_sq  ){
  weights <- k_j/ (k_j*tau_sq + (k_j-1)*rho*sigma_j_sq + omega_sq + sigma_j_sq)
  
  return(weights)
  
  
}


# degrees of freedom for moderator with multiple categories
multiple_categories <- function(data= NULL, 
                                moderator_val, 
                                cluster_id, 
                                sigma_j_sq_val, 
                                rho_val, 
                                omega_sq_val, 
                                tau_sq_val){
  
  ### get variables
  
  if (!is.null(data)) {
    moderator_call <- substitute(moderator_val)
    cluster_call <- substitute(cluster_id)
    sigma_j_sq_call <- substitute(sigma_j_sq_val)
    rho_call <- substitute(rho_val)
    omega_sq_call <- substitute(omega_sq_val)
    tau_sq_call <- substitute(tau_sq_val)
    
    
    env <- list2env(data, parent = parent.frame())
    
    moderator_val <- eval(moderator_call, env)
    cluster_id <- eval(cluster_call, env)
    sigma_j_sq_val <- eval(sigma_j_sq_call, env)
    rho_val <- eval(rho_call, env)
    omega_sq_val <- eval(omega_sq_call, env)
    tau_sq_val <- eval(tau_sq_call, env)
    
  }
  
  dat <- data.frame(
    cluster = factor(cluster_id),
    moderator = moderator_val,
    sigma_j_sq = sigma_j_sq_val,
    sigma_j_sq = sigma_j_sq_val,
    rho = rho_val,
    omega_sq = omega_sq_val, 
    tau_sq = tau_sq_val
    
    
    
  )
  
  
  c <-  length(unique(dat$moderator))
  q <-  c-1
  
  df_num <- q
  
  study_level_data <-  dat %>% 
    group_by(cluster) %>% 
    mutate(k_j = n()) %>% ungroup()

  
  
  by_category <-  study_level_data %>% 
    group_by(moderator) %>% 
    mutate(
      w_j = study_level_weights(k_j =k_j ,sigma_j_sq = sigma_j_sq, tau_sq = tau_sq, rho = rho, omega_sq = omega_sq )) %>% 
    summarise(W = sum(w_j),
              E_Vr = 1/W,
              nu_Vr = (sum(w_j^2/(W-w_j)^2) - ((2/W)*sum(w_j^3/(W-w_j)^2)) + ((1/W^2)*(sum(w_j^2/(W-w_j)))^2)) ^(-1), .groups = 'drop' )
  
  
  
  by_category$nu_D <- (q*(q+1)) / (2*(sum( (1/by_category$nu_Vr)* (1-(1/(by_category$E_Vr*sum(by_category$W))))^2 )))
  by_category$df_num = df_num
  
  
  return(by_category)
  
}

#NCP
ncp <- function(weights, mu_p){
  
  mu_wavg = sum(weights*(mu_p))/ sum(weights)
  ncp = sum(weights*(mu_p-mu_wavg)^2)
  return(ncp)
}


# power for the CHE-RVE model

power_CHE_RVE_study_cat <- function(data = NULL, 
                                    moderator_val, 
                                    cluster_id, 
                                    sigma_j_sq_val, 
                                    rho_val, 
                                    omega_sq_val, 
                                    tau_sq_val,
                                    mu, 
                                    alpha) {
  
  ### get variables
  
  if (!is.null(data)) {
    moderator_call <- substitute(moderator_val)
    cluster_call <- substitute(cluster_id)
    sigma_j_sq_call <- substitute(sigma_j_sq_val)
    rho_call <- substitute(rho_val)
    omega_sq_call <- substitute(omega_sq_val)
    tau_sq_call <- substitute(tau_sq_val)
    
    
    env <- list2env(data, parent = parent.frame())
    
    moderator_val <- eval(moderator_call, env)
    cluster_id <- eval(cluster_call, env)
    sigma_j_sq_val <- eval(sigma_j_sq_call, env)
    rho_val <- eval(rho_call, env)
    omega_sq_val <- eval(omega_sq_call, env)
    tau_sq_val <- eval(tau_sq_call, env)
    
  }
  
  dat <- data.frame(
    cluster = factor(cluster_id),
    moderator = moderator_val,
    sigma_j_sq = sigma_j_sq_val,
    sigma_j_sq = sigma_j_sq_val,
    rho = rho_val,
    omega_sq = omega_sq_val, 
    tau_sq = tau_sq_val
    
    
    
  )
  
  
  
  df <- multiple_categories(dat = dat, moderator = moderator, cluster = cluster,sigma_j_sq = sigma_j_sq, rho = rho,omega_sq = omega_sq, tau_sq = tau_sq)
  
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

#code for cdf of non-central f dist.
# non_central_f_cdf_j <- function( lambda, J,  d1, d2, x){
#   
#   f  <- (((((1/2)*lambda)^J)/factorial(J))*(exp(1))^(-lambda/2)*pbeta(((d1*x)/(d2+d1*x)), (d1/2 + J), (d2/2))) 
#   
# }
# 
# 
# non_central_f_cdf <- function(J, lambda,  d1, d2, x){
#   
#   f <-sapply(J, non_central_f_cdf_j,lambda = lambda, d1 = d1, d2 = d2, x = x)
#   sum(unlist(f))
#   
# }
# 
# 
# non_central_f_cdf_reverse <- function(J, lambda,  d1, d2, x, area){
#   
#   f <-sapply(J, non_central_f_cdf_j,lambda = lambda, d1 = d1, d2 = d2, x = x)
#   sum(unlist(f)) - area
#   
# }

# find lambda (NCP)

find_lambda <- function(lambda, d1, d2, x, area, interval, tol){
  
 roots <- uniroot(f = non_central_f_cdf_reverse, d1 = d1, d2 = d2, x = x, area = area, interval = interval,  extendInt = "yes", tol = tol)
  
 return(roots$root)
 
}


# patterns of the beta_coefficients
f_c <- function(pattern){
  
  # C  = 2
  if(pattern == 1){
    f_c <- c(0, 1)
  }
  
  # C  = 3
  if(pattern == 2){
    f_c <- c(0, 1, 2)
  }
  
  if(pattern == 3){
    f_c <- c(0, 0, 1)
  }
  
  if(pattern == 4){
    f_c <- c(0, 1, 1)
  }
  
  # C  = 4
  if(pattern == 5){
    f_c <- c(0, 1, 2, 3)
  }
  
  if(pattern == 6){
    f_c <- c(0, 0, 1, 2)
  }
  
  if(pattern == 7){
    f_c <- c(0, 0, 0, 1)
  }
  
  if(pattern == 8){
    f_c <- c(0, 0, 1, 1)
  }
  
  return(f_c)
}







# find the scalling factor
zeta <- function(pattern, lambda, weights){
  
  zeta <- sqrt((lambda)/(sum(weights*(pattern - sum(weights*pattern)/sum(weights) )^2)))
  
  return(zeta)
}

# beta coefficients 
build_mu <- function(pattern, zeta){
  
  pattern*zeta
  
}

# design matrix 
design_matrix <- function(C, J, bal){

  
  if(bal == "balanced"){
    
    J_c <- J/C
    
  }
  
  
  categories <-   data.frame(cat = rep(LETTERS[1:C], each = J_c))
  
  X <- model.matrix(~ 0 + cat, categories) 
  
  return(X)
  
}

#X <- design_matrix(C= 4, J= 24, bal = "balanced")


# data for approximation function tests
dat_approx <- function(C, J, tau_sq, omega_sq, rho, k_j, n_j) {
  
  J_c <- J/C
 
  
  dat_approx <- tibble(studyid = c(1:J),
                       k_j = rep(k_j, J),
                       n_j = rep(n_j, J),
                       sigma_j_sq = 4 / n_j, #### double check this
                       omega_sq = rep(omega_sq, J),
                       rho = rep(rho, J),
                       tau_sq = rep(tau_sq, J),
                       cat = rep(LETTERS[1:C], each = J_c)
  )
  
  return(dat_approx)
}

