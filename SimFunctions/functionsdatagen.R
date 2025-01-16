

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
  
  study_level_data <-  dat |> 
    group_by(cluster) |>
    mutate(k_j = n()) |> ungroup()

  
  
  by_category <-  study_level_data |> 
    group_by(moderator) |> 
    mutate(
      w_j = study_level_weights(k_j =k_j ,sigma_j_sq = sigma_j_sq, tau_sq = tau_sq, rho = rho, omega_sq = omega_sq )) |> 
    summarise(W = sum(w_j),
              E_Vr = 1/W,
              nu_Vr = (sum(w_j^2/(W-w_j)^2) - ((2/W)*sum(w_j^3/(W-w_j)^2)) + ((1/W^2)*(sum(w_j^2/(W-w_j)))^2)) ^(-1), .groups = 'drop' )
  
  
  
  by_category$nu_D <- (q*(q+1)) / (2*(sum( (1/by_category$nu_Vr)* (1-(1/(by_category$E_Vr*sum(by_category$W))))^2 )))
  by_category$df_num = df_num
  
  
  return(by_category)
  
}

#  NCP
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







# find the scaling factor
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

mu_values <- function(
    #C, 
  J, 
  tau_sq, 
  omega_sq, 
  rho, 
  P, 
  k_j = 3, ####CHANGE THIS LATER
  n_j = 30, ####CHANGE THIS LATER
  f_c_val ) {
  
  # C  = 2
  if(f_c_val == 1){
    C= 2
  }
  
  # C  = 3
  if(f_c_val == 2){
    C= 3
  }
  
  if(f_c_val == 3){
    C= 3
  }
  
  if(f_c_val == 4){
    C= 3
  }
  
  # C  = 4
  if(f_c_val == 5){
    C= 4
  }
  
  if(f_c_val == 6){
    C= 4
  }
  
  if(f_c_val == 7){
    C= 4
  }
  
  if(f_c_val == 8){
    C= 4
  }
  
  
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
  
  
  possibly_find_lambda <- possibly(.f = find_lambda, otherwise = NA)
  
  
  lambda <- possibly_find_lambda(d1= df_num, 
                                 d2 = df_den, 
                                 x = qf(1-.05, df_num,df_den ), 
                                 area = 1 - P, 
                                 interval = c(0, 100000), 
                                 tol = 0.0001)
  
  if(!is.na(lambda)) {
    zeta_val <- zeta(pattern =  f_c(pattern = f_c_val), lambda = lambda, weights =dfs$W)
    
    mu_values <- build_mu(pattern =  f_c(pattern = f_c_val), zeta = zeta_val)
    
  }else{
    mu_values <- NA
  }
  
  
  
  return(mu_values)
  
  
}

#-----------------------------------------------------------------------------
### functions below are pulled and modified from https://osf.io/gaz9t/ (the study
### I'm extending to be consistent in workflow. 
### probably try to make variable names consistent at some point...
#-----------------------------------------------------------------------------

#smd -- pulled and slightly modified from https://osf.io/e4npa

generate_smd <- function(delta, k, N, Sigma) {
  
  # make sure delta is a vector
  delta_vec <- rep(delta, length.out = k)
  
  # create Sigma matrix assuming equicorrelation
  if (!is.matrix(Sigma)) Sigma <- Sigma + diag(1 - Sigma, nrow = k) # cor matrix for 1 study
  
  # generate numerator of SMD 
  mean_diff <- rmvnorm(n = 1, mean = delta_vec, sigma = (4/N) * Sigma) 
  
  # covariance 
  cov_mat <- as.matrix(rWishart(n = 1, df = N - 2, Sigma = Sigma)[,,1])
  sigma_sq <- diag(cov_mat) / (N - 2)
  
  # SMD
  d <- as.vector(mean_diff / sqrt(sigma_sq))  # cohen's d 
  bias_corr <- (1 - (3/((4 * (N - 2)) - 1)))
  g <- d * bias_corr # Hedges g
  var_g <- bias_corr^2 * (4 / N + d^2 / (2 * (N - 2)))
  
  dat <- tibble(g = g, var_g = var_g)
  
  return(dat)
}

#-----------------------------------------------------------------------------  
### we are fixing rho to be the same across studies, but I'll keep this 
### function here in case that changes. they obtained this reparameterization 
### function from https://osf.io/gaz9t/ 

# reparm <- function(cor_sd, cor_mu) {
#   alpha <- cor_mu * ((cor_mu * (1-cor_mu) / cor_sd^2 )- 1)
#   bet <-  (1-cor_mu) * ((cor_mu *(1-cor_mu) / cor_sd^2 )- 1)
#   reparms <- as.data.frame(cbind(alpha, bet))
#   return(reparms)
# }
#-----------------------------------------------------------------------------  

generate_meta <- function(J, tau_sq, 
                          omega_sq, bal, C,
                         # cor_mu, cor_sd, 
                          rho, P, sample_sizes, k_j, 
                         ### added k_j here for now, but it will probably
                         # be part of sample_sizes so should remove later. 
                          nj, f_c_val,
                          return_study_params = FALSE,
                          seed = NULL) {
  
  require(dplyr)
  require(purrr)
  if (!is.null(seed)) set.seed(seed)
  
  # Study data --------------------------------------------------------------
  
  # cor_params <- reparm(cor_sd=cor_sd, cor_mu=cor_mu)
  mu_vector <- mu_values(J = J, tau_sq = tau_sq, omega_sq = omega_sq, 
                         rho = rho, P = P, k_j = k_j, n_j = nj, f_c_val = f_c_val) 
  #####3may need to replace k_j input with sample_sizes later..
  
  X <- design_matrix(C= C, J= J, bal = bal)
  
  Xbeta <- as.numeric(X %*% mu_vector)
  
  study_data <- 
    tibble(
      N = rep(sample_sizes, length.out = J),
      nj = rep(nj, length.out = J), 
    #  Sigma = rbeta(n=J, shape1=cor_params$alpha, shape2=cor_params$bet),
      Sigma = rep(rho, J) ## assume equi-correlation across studies -- 
      #############ask James about this one
    ) %>%
    mutate(
      u_j = rnorm(J, 0, sqrt(tau_sq)),
      v_ij = map(nj, ~ rnorm(., 0, sqrt(omega_sq))),
      delta = map2(u_j, v_ij, ~ Xbeta + .x + .y)
     # delta2 = Xbeta + u_j + v_ij, 
    ) %>%
    select(delta, nj, N, Sigma)
  
  if (return_study_params) return(study_data)
  
  # Generate full meta data  -----------------------------------------------
  
  # first line runs generate_smd
  meta_reg_dat <- 
    pmap_df(study_data, generate_smd, .id = "studyid") %>%
    mutate(
      esid = 1:n()
    )
  
  
  return(meta_reg_dat)
}

# Function for random sampling n rows from dataset

n_ES_empirical <- function(dat, J, with_replacement = TRUE) {
  dat[sample(NROW(dat), size = J, replace = with_replacement),]
}

