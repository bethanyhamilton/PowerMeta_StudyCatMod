

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
    rho = rho_val,
    omega_sq = omega_sq_val, 
    tau_sq = tau_sq_val
  )
  
  c <-  length(unique(dat$moderator))
  q <-  c - 1
  
  study_level_data <-  dat |> 
    group_by(moderator) |>
    mutate(mean_var = mean(sigma_j_sq), k_j = n()) |> 
    ungroup() |> 
    select(-sigma_j_sq)
  
  
  by_category <-  
    study_level_data |> 
    group_by(moderator) |> 
    mutate(
      w_j = study_level_weights(
        k_j = k_j,
        sigma_j_sq = mean_var, 
        tau_sq = tau_sq, 
        rho = rho, 
        omega_sq = omega_sq 
      )
    ) |> 
    summarise(
      W = sum(w_j),
      E_Vr = 1/W,
      nu_Vr = (sum(w_j^2 / (W - w_j)^2) - (2 / W) * sum(w_j^3 / (W - w_j)^2) + (1 / W^2) * sum(w_j^2 / (W - w_j))^2)^(-1), 
      .groups = 'drop' 
    )
  
  
  return(by_category)
  
}


# function for just df:

df_val <- function(E_Vr, nu_Vr, W) { 
  
  q <- length(E_Vr) - 1
  
  df_num <-  q 
  nu_D <- q * (q+1)  / (2 * sum( (1 / nu_Vr) * (1 - (1 / (E_Vr * sum(W))))^2 ))
   
  
  return(tibble(df_num, nu_D))
  
  }
   



# non-centrality parameter -- NCP
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
  

  #rename functions and objects 
  multiple_cat <- multiple_categories(
                          data = data, 
                          moderator_val = moderator_val, 
                          cluster_id = cluster_id, 
                          sigma_j_sq_val = sigma_j_sq_val, 
                          rho_val = rho_val,
                          omega_sq_val = omega_sq_val, 
                          tau_sq_val = tau_sq_val
                          )
  

  
  
  ncp <-  ncp(weights = multiple_cat$W, mu_p = mu)
  
  
  ##########################################
  
  
  
  df_num <- df_val(E_Vr = multiple_cat$E_Vr, nu_Vr = multiple_cat$nu_Vr, W = multiple_cat$W)[[1]]
    
    
  #  df$df_num[1]
  df_den <- df_val(E_Vr = multiple_cat$E_Vr, nu_Vr = multiple_cat$nu_Vr, W = multiple_cat$W)[[2]] - df_val(E_Vr = multiple_cat$E_Vr, nu_Vr = multiple_cat$nu_Vr, W = multiple_cat$W)[[1]] + 1
  
  
  F_crit <- qf(1-alpha, df_num, df_den)
  power <- 1- pf(F_crit, df_num, df_den, ncp = ncp)
  
  
  power_ncp <- data.frame(power, ncp, df_den)
  
  return(power_ncp)
  
  
}





# plug into uniroot to get lambda or the non-centrality parameter
non_central_f_cdf_reverse <- function(lambda, d1, d2, x, area) {
  
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

find_lambda <- function(d1, d2, x, area, interval, tol){
  
  
  # needed to add this since mu should be 0 in this case since it is equal to
  # type I error of .05 (we are only looking at that condition)
  if (area == 0.95) {
    return(0)
  }
  
 roots <- uniroot(f = non_central_f_cdf_reverse, d1 = d1, 
                  d2 = d2, x = x, area = area, interval = interval,  
                  extendInt = "yes", tol = tol)
  
 return(roots$root)
 
}

# patterns of the beta_coefficients
f_c_lookup <- list(
  P1 = c(0,1),
  P2 = c(0,1,2),
  P3 = c(0, 0, 1),
  P4 = c(0, 1, 1),
  P5 = c(0, 1, 2, 3),
  P6 = c(0, 0, 1, 2),
  P7 = c(0, 0, 0, 1),
  P8 = c(0, 0, 1, 1)
)


# find the scaling factor
zeta <- function(pattern, lambda, weights){
  
  zeta <- sqrt(lambda / sum(weights*(pattern - sum(weights * pattern) / sum(weights))^2))
  
  return(zeta)
}


# design matrix 

### options for balance of the number of studies across categories:
#### 1. set scenarios where each J_c add up a multiple of J = 12
#### 2. fractions in combination with multiple of J = 12 that result in whole numbers
#### 3. stochastic assignment  --  where each category has a probability of getting assigned to a study

#### will go with 2 since we can see what is going on exactly with different combinations of factors. 
#### constrained to multiples of 12.to compare across balanced and unbalanced conditions -- probably need to workshop this more
#follow rule -- always increasing in number of studies. smaller number of studies will be on smallest mu. 

design_matrix <- function(C, J, bal, k_j){

  categories <- mod(C = C, J = J, bal = bal, k_j = k_j)
  
  X <- model.matrix(~ 0 + category, categories) 
  
  return(X)
  
}

# data set with categorical moderator
mod <- function(C, J, bal, k_j ){

  N = sum(k_j)
  
  
 category_study  <-  dat_approx(C = C, 
                                J = J, 
                                tau_sq = NULL, 
                                omega_sq = NULL, 
                                rho = NULL, 
                                k_j = k_j, 
                                N= NULL, 
                                sigma_j_sq = NULL, 
                                bal = bal)
  
  
  #needs to be K x K
  
  covariates <- tibble(category = c(rep(category_study$cat,  k_j)),
         studyid = as.character(c(rep(c(1:J),  k_j))),
         esid = 1:N)

  
  return(covariates)
  
}

#study-level dat
dat_approx <- function(C, J, tau_sq, omega_sq, rho, k_j, N= NULL, sigma_j_sq = NULL, bal) {
  
  

  if(!is.null(N) & is.null(sigma_j_sq)){
    sigma_j_sq = 4 / N
    
    
  }
    
  if(bal == "balanced_j"){
    
    J_c <- J/C
    cat <-   rep(LETTERS[1:C], each = J_c)
    
    
  }
  
  ##### unbalanced
  if(bal == "unbalanced_j"){
    
    if (C == 2){
      J_vec <- J * c(1/4, 3/4)
      cat <-   rep(LETTERS[1:C],  J_vec)
    }
    
    # #less unbalanced
    # if (C == 3){
    #   
    #   J_vec <- J * c(1/6, 5/12, 5/12)
    #   
    #   
    # }
    
    if (C == 3){
      J_vec <- J * c(1/4, 1/4, 1/2)
      cat <-   rep(LETTERS[1:C],  J_vec)
    }
    
    
    if (C == 4){
      J_vec <- J * c(1/6, 1/6, 1/6, 1/2)
      cat <-   rep(LETTERS[1:C],  J_vec)
    }
  }
 

dat_approx <- tibble(studyid = c(1:J),
                       k_j= k_j,
                       sigma_j_sq = sigma_j_sq, ##### CHANGE THIS LATER
                       omega_sq = rep(omega_sq, J),
                       rho = rep(rho, J),
                       tau_sq = rep(tau_sq, J),
                       cat = cat
  )
  
  return(dat_approx)
}




# function to get power value given a set of beta coefficients/mu values
power_approximation <- function(
     
                                J, 
                                tau_sq, 
                                omega_sq, 
                                rho, 
                                N_mean = NULL,
                                k_mean = NULL,
                                sigma_j_sq_mean = NULL,
                                N_dist = NULL,
                                sigma_j_sq_dist = NULL,
                                pilot_data = NULL,
                                iterations = 1L,
                                sample_size_method = c("balanced","stylized","empirical"),
                                P, 
                                f_c_val,  
                                bal,
                                seed = NULL
                                
){
  
  
  
  
  # -----------------------------------------------------
  ## sampling methods code below pulled from https://osf.io/gaz9t/  and 
  ## incorporated into this study
  
  ### add in sample sigma_j_sq*************************
  # -------------------------------------------------------  
  
  if (!is.null(seed)) set.seed(seed)
  N <- list()
  kjs <- list()
  sigma_j_sqs <- list()
  
  # Average sample sizes and number of effect sizes
  if ("balanced" %in% sample_size_method) {
    
    if (is.null(N_mean) | is.null(k_mean)) stop("Must specify values for N_mean and k_mean.")
    
    N <- c(N, balanced = N_mean)
    kjs <- c(kjs, balanced = k_mean)
    sigma_j_sqs <- c(sigma_j_sqs, balanced = sigma_j_sq_mean)
  }
  
  # Stylized sample size and n ES distributions
  if ("stylized" %in% sample_size_method) {
    
    if (is.null(N_dist) | is.null(k_mean)) stop("Must specify values for sigma2_dist and k_mean.")
    
    styled_Ns <- map(1:iterations, ~pmax(10, rgamma(J, shape = N_dist$estimate[1], rate = N_dist$estimate[2])))
    styled_sigma_j_sqs <- map(1:iterations, ~pmax(0.008, rgamma(J, shape = sigma_j_sq_dist$estimate[1], rate = sigma_j_sq_dist$estimate[2])))
    styled_kjs <- map(1:iterations, ~ 1 + rpois(J, k_mean - 1))
    
    N <- c(N, stylized = styled_Ns)
    kjs <- c(kjs, stylized = styled_kjs)
    sigma_j_sqs <- c(sigma_j_sqs, balanced = styled_sigma_j_sqs)
    
  }
  
  
  # Empirical sample sizes and k ES distributions
  if ("empirical" %in% sample_size_method) {
    
    if (is.null(pilot_data)) stop("Must specify a dataset with pilot_data.")
    
    pilot_sample <- map(1:iterations, ~ n_ES_empirical(pilot_data, J))
    N <- c(N, empirical = map(pilot_sample, ~ .x$N))
    kjs <- c(kjs, empirical =  map(pilot_sample, ~ .x$kj))
    sigma_j_sqs <- c(sigma_j_sqs, empirical =  map(pilot_sample, ~ .x$sigma_j_sq))
  }
  
  
  C <-  ftoc[[f_c_val]]
  

  # res_CHE_RVE <- map2_dfr(
  #   .x = N, .y = kjs, .f = run_power, C=C,
  #   J = J, tau_sq = tau_sq, omega_sq = omega_sq,  rho = rho, bal = bal, P =P, f_c_val = f_c_val,
  #   .id = "samp_method"
  # )
  
  res_CHE_RVE <- pmap_df(
    list( kjs, N, sigma_j_sqs), .f = run_power, C=C,
       J = J, tau_sq = tau_sq, omega_sq = omega_sq,  rho = rho, bal = bal, P =P, f_c_val = f_c_val,
       .id = "samp_method"
  )
  
  return(res_CHE_RVE)
  
}


# JEP: Do you actually need this function?
# It seems like it does basically the same thing as power_CHE_RVE_study_cat()?
run_power <- function(C, 
                      J, 
                      tau_sq, 
                      omega_sq, 
                      rho, 
                      k_j, 
                      P, 
                      f_c_val,
                      bal,
                      N = NULL, 
                      sigma_j_sq = NULL
){
  

  mu_vec <- mu_values(J = J, tau_sq = tau_sq, 
                      omega_sq = omega_sq,
                      rho = rho, P = P, 
                      k_j = k_j, N = N, 
                      sigma_j_sq = sigma_j_sq, 
                      f_c_val = f_c_val, 
                      bal = bal)
  
  dat_app <-  dat_approx(C = C, J = J, tau_sq = tau_sq, 
                         omega_sq = omega_sq, rho = rho, 
                         k_j = k_j, N = N , sigma_j_sq = sigma_j_sq ,bal = bal)
  
  
  power <-  power_CHE_RVE_study_cat(data = dat_app, 
                                    moderator_val  = dat_app$cat, 
                                    cluster_id = dat_app$studyid, 
                                    sigma_j_sq_val = dat_app$sigma_j_sq,
                                    rho_val = dat_app$rho,
                                    omega_sq_val = dat_app$omega_sq,
                                    tau_sq_val = dat_app$tau_sq,
                                    mu = mu_vec, 
                                    alpha = .05)
  
  
  
  return(power)
  
}

# list of pattern to number of categories C
ftoc <- c(2, 3, 3, 3, 4, 4, 4, 4)
names(ftoc) <- c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8")


# given conditions get beta coefficients/ mu values
mu_values <- function(
  J, 
  tau_sq, 
  omega_sq, 
  rho, 
  P, 
  k_j,
  N = NULL, 
  sigma_j_sq = NULL,
  bal,
  f_c_val ) {
  

  
  C <- ftoc[[f_c_val]]

  
  dat_app <-  dat_approx(C = C, J = J, tau_sq = tau_sq, 
                         omega_sq = omega_sq, rho = rho, 
                         k_j = k_j,
                         N = N, 
                         sigma_j_sq = sigma_j_sq,
                         bal = bal
                         )
  
  
  
  multiple_cat <- multiple_categories(
    
    data = dat_app,
    moderator_val = cat, 
    cluster_id = studyid, 
    sigma_j_sq_val = sigma_j_sq, 
    rho_val = rho,
    omega_sq_val = omega_sq, 
    tau_sq_val = tau_sq
  )
  
  df_num <- df_val(E_Vr = multiple_cat$E_Vr, nu_Vr = multiple_cat$nu_Vr, W = multiple_cat$W)[[1]]
  df_den <- df_val(E_Vr = multiple_cat$E_Vr, nu_Vr = multiple_cat$nu_Vr, W = multiple_cat$W)[[2]] - df_val(E_Vr = multiple_cat$E_Vr, nu_Vr = multiple_cat$nu_Vr, W = multiple_cat$W)[[1]] + 1
  
  
  
  possibly_find_lambda <- possibly(.f = find_lambda, otherwise = NA)
  
  
  lambda <- possibly_find_lambda(d1= df_num, 
                                 d2 = df_den, 
                                 x = qf(1-.05, df_num,df_den ), 
                                 area = 1 - P, 
                                 interval = c(0, 100000), 
                                 tol = 0.0001)
  
  if(!is.na(lambda)) {
    zeta_val <- zeta(
                     pattern = f_c_lookup[[f_c_val]], 
                     lambda = lambda, weights =multiple_cat$W)
    
    # beta coefficients 
    mu_values <-  f_c_lookup[[f_c_val]] * zeta_val
      
    
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
# Function for random sampling n rows from dataset
n_ES_empirical <- function(dat, J, with_replacement = TRUE) {
  dat[sample(NROW(dat), size = J, replace = with_replacement),]
}


#smd -- pulled and slightly modified from https://osf.io/e4npa
generate_smd <- function(delta, k_j, N, Sigma) {
  
  # make sure delta is a vector
  delta_vec <- rep(delta, length.out = k_j)
  
  # create Sigma matrix assuming equicorrelation
  if (!is.matrix(Sigma)) Sigma <- Sigma + diag(1 - Sigma, nrow = k_j) # cor matrix for 1 study
  
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

# generate_meta <- function(return_mu = NULL, J, tau_sq, 
#                           omega_sq, bal, 
#                           rho, P, sample_sizes, k_j,
#                           f_c_val,
#                           sigma_j_sq = NULL,
#                           return_study_params = FALSE,
#                           seed = NULL) {
#   
#   require(dplyr)
#   require(purrr)
#   if (!is.null(seed)) set.seed(seed)
#   
#   # Study data --------------------------------------------------------------
#   
#   C <-  ftoc[[f_c_val]]
#   
#   # cor_params <- reparm(cor_sd=cor_sd, cor_mu=cor_mu)
#   mu_vector <- mu_values(J = J, tau_sq = tau_sq, omega_sq = omega_sq, 
#                          rho = rho, P = P, k_j = k_j,  f_c_val = f_c_val,
#                          bal =bal, sigma_j_sq = sigma_j_sq, N = sample_sizes) 
#   
#   X <- design_matrix(C= C, J= J, bal = bal, k_j = k_j)
#   
#   Xbeta <- as.numeric(X %*% mu_vector)
#   
#   mod_data = mod(C= C, J= J, bal = bal, k_j = k_j) |> 
#     mutate(studyid = factor(x = studyid, levels = 1:J))
#   
#   mod_data_stud <- mod_data |> 
#     group_nest(studyid, .key = "X")
#   
#   n_ES_total <- nrow(X)
#   
#   
#   
#   u_j = rnorm(J, 0, sqrt(tau_sq))[mod_data$studyid]
#   v_ij = rnorm(n_ES_total, 0, sqrt(omega_sq))
#   
#   study_data <- 
#     tibble(
#       delta = Xbeta + u_j + v_ij,
#       studyid = mod_data$studyid
#     ) |>
#     group_by(studyid) |>
#     summarize(
#       delta = list(delta),
#       k_j = n(),
#       .groups = "drop"
#     ) |>
#     mutate(
#       #  Sigma = rbeta(n=J, shape1=cor_params$alpha, shape2=cor_params$bet),
#       Sigma = rep(rho, J),
#       N = sample_sizes,
#     ) 
#   
#   
#     
#   if (return_study_params) return(study_data)
#   
#   # Generate full meta data  -----------------------------------------------
#   
#   meta_reg_dat <- 
#     study_data |> 
#     {\(.) dplyr::mutate(.,
#       smds = pmap(select(., -studyid), generate_smd)
#     )}() |>
#     left_join(mod_data_stud, by = "studyid") |>
#     select(-delta, -k_j, -N, -Sigma) |>
#     unnest(cols = c(smds, X))
#   
#   if(!is.null(return_mu)){
#     meta_reg_dat <- list(meta_reg_dat, mu_vector)
#   }
#   
#   
#   return(meta_reg_dat)
# }



generate_meta <- function(J, tau_sq, 
                          omega_sq, bal, 
                          mu_vector,
                          rho, P, 
                          sample_sizes, k_j,
                          f_c_val,
                          sigma_j_sq = NULL,
                          return_study_params = FALSE,
                          seed = NULL) {
  
  require(dplyr)
  require(purrr)
  if (!is.null(seed)) set.seed(seed)
  
  # Study data --------------------------------------------------------------
  
  C <-  ftoc[[f_c_val]]
  
  # cor_params <- reparm(cor_sd=cor_sd, cor_mu=cor_mu)
  
  X <- design_matrix(C= C, J= J, bal = bal, k_j = k_j)
  
  Xbeta <- as.numeric(X %*% mu_vector)
  
  mod_data = mod(C= C, J= J, bal = bal, k_j = k_j) |> 
    mutate(studyid = factor(x = studyid, levels = 1:J))
  
  mod_data_stud <- mod_data |> 
    group_nest(studyid, .key = "X")
  
  n_ES_total <- nrow(X)
  
  
  
  u_j = rnorm(J, 0, sqrt(tau_sq))[mod_data$studyid]
  v_ij = rnorm(n_ES_total, 0, sqrt(omega_sq))
  
  study_data <- 
    tibble(
      delta = Xbeta + u_j + v_ij,
      studyid = mod_data$studyid
    ) |>
    group_by(studyid) |>
    summarize(
      delta = list(delta),
      k_j = n(),
      .groups = "drop"
    ) |>
    mutate(
      #  Sigma = rbeta(n=J, shape1=cor_params$alpha, shape2=cor_params$bet),
      Sigma = rep(rho, J),
      N = sample_sizes,
    ) 
  
  
  
  
  if (return_study_params) return(study_data)
  
  # Generate full meta data  -----------------------------------------------
  
  meta_reg_dat <- 
    study_data |> 
    {\(.) dplyr::mutate(.,
                        smds = pmap(select(., -studyid), generate_smd)
    )}() |>
    left_join(mod_data_stud, by = "studyid") |>
    select(-delta, -k_j, -N, -Sigma) |>
    unnest(cols = c(smds, X))
  
  
  
  
  return(meta_reg_dat)
}



# --------------------------------------------------------------------------
## Estimate Model -- straightforward we are not looking at mispecification at
## the moment. Just need to run CHE-RVE model
# --------------------------------------------------------------------------




estimate_model <- function(data = NULL,
                          # return_mu = NULL,
                           moderator_val,
                           cluster_id,
                           delta, 
                           delta_var,
                           es_id,
                           r= 0.7,
                           smooth_vi = TRUE, 
                           control_list = list()
){
  
  require(dplyr)
  
  
  # if(!is.null(return_mu)){
  #  
  #    dat <- data[[1]]
  #    mu_vector <- data[[2]]
  #    
  #    
  # }
    
  dat <- data.frame(
    study_id = cluster_id,
    moderator = moderator_val,
    g = delta,
    vi = delta_var,
    esid = es_id
    
    
    
  )
  

  if (smooth_vi) { 
    dat <- 
      dat |> 
      group_by(study_id) |>
      mutate(var_g_j = mean(vi, na.rm = TRUE)) |>
      ungroup()
  }
 
  C = length(unique(dat$moderator))
  
  
  V_list <- 
    vcalc(
      vi = var_g_j,
      cluster = study_id,
      rho = r,
      obs = esid,
      sparse = TRUE,
      data = dat 
    )
  
  rma_fit <- 
    purrr::possibly(metafor::rma.mv, otherwise = NULL)(
      g ~ 0 + moderator,
      V = V_list, 
      random = ~ 1 | study_id / esid,
      data = dat,
      test = "t",
      sparse = TRUE,
      verbose = FALSE,
      control = control_list 
    )
  
  tryCatch({ 
  coef_RVE <-  robust(
    rma_fit, # estimation model above
    cluster = study_id, # define clusters
    clubSandwich = TRUE # use CR2 adjustment
  )
  
  wald_test_results <- Wald_test((rma_fit), constraints = constrain_equal(1:C), vcov =  "CR2")
  
  
  
  res <- 
    tibble(
      est = list(coef_RVE$beta),
      est_var = list(coef_RVE$se^2),
      df1 = wald_test_results$df_num,
      df2 = wald_test_results$df_denom,
      p_val = wald_test_results$p_val, 
      model = "CHE",
      var = "RVE"
    )
  
 
 
  

  }, error = function(w) { res<- 
    tibble(
      est = NA_real_,
      est_var = NA_real_,
      df1 = NA_real_,
      df2 = NA_real_,
      p_val = NA_real_, 
      model = "CHE",
      var = "RVE"
    )
  
  })
  
  # if(!is.null(return_mu)){
  #   
  #   res$mu_vector_list <- list(mu_vector)
  #   
  #   
  # } 
  res
  
}



#-----------------------------------------------------------
#### Simulation performance
#-----------------------------------------------------------

#### add in est_M, est_V, var_M, var_df for the mu_estimates 





# sim_performance <- function(results) {
# 
#   require(dplyr)
# 
#   results %>%
#     summarise(
#       n_sim = n(),
#       cnvg = mean(!is.na(p_val)),
#     #  df1_mean = mean(df1),
#     #  df2_mean = mean(df2),
#       rej_rate_05 = mean(p_val < 0.05, na.rm = TRUE)
#     )
# 
# }


# Approximate vs. simulated power
#### add function that joins the approximation results and sim results and 
#### find the deviance between the two


#-----------------------------------------------------------
#### Simulation Driver
#-----------------------------------------------------------

## old one with the power conditional on study features (mu generated in generate_meta())
# run_sim <- function(iterations,
#                     J, 
#                     tau_sq, 
#                     omega_sq, 
#                     bal, 
#                     rho, 
#                     P, 
#                     f_c_val,
#                     sigma_j_sq_inc = NULL,
#                     pilot_data = NULL,
#                     return_study_params = FALSE,
#                     seed = NULL,
#                     summarize_results = FALSE){
#   
#  require(dplyr)
#   
#   #in this case power is conditional on study features. basically I specify given set of study features that are
#   #generated randomly, and I want the prob of sig. result to be a value for this given set of study features. 
#   #controlling power given  x means if you specify a different x the mu values will need to change to hold power constant
#   #which is a little perplexing
# 
#  if (!is.null(seed)) set.seed(seed)
#   
# 
#  
#  results <- map(1:iterations, ~{
#    
#                   sample_dat <- n_ES_empirical(pilot_data, J = J)
#                   
#                   if(sigma_j_sq_inc == TRUE){
#                     sigma_j_sq = sample_dat$sigma_j_sq
#                   } else{
#                     sigma_j_sq = NULL
#                   }
#                   
#                   dat <- generate_meta(J = J, 
#                                        tau_sq = tau_sq, 
#                                        omega_sq = omega_sq, 
#                                        bal = bal, 
#                                        rho = rho, 
#                                        P = P, 
#                                        sample_sizes = sample_dat$N, 
#                                        k_j = sample_dat$kj,
#                                        f_c_val = f_c_val,
#                                        sigma_j_sq = sigma_j_sq,
#                                        return_study_params = return_study_params,
#                                        return_mu = NULL,
#                                        seed = seed ## maybe should remove this here and only keep in run_sim
#                                      )
#                   
#               
#                   est_res <-  estimate_model(data = dat,
#                                              moderator_val = dat$category,
#                                              cluster_id = dat$studyid,
#                                              delta = dat$g, 
#                                              delta_var =  dat$var_g,
#                                              es_id = dat$esid,
#                                              r= rho,
#                                              smooth_vi = TRUE, 
#                                              return_mu = NULL,
#                                              control_list = list()
#                   )
#                   
#                   }) |> dplyr::bind_rows()
#   
#  if (summarize_results) {
#    performance <- sim_performance(results = results)
#    return(performance)
#  } else {
#    return(results)
#  }
# }

run_sim <- function(iterations,
                    J, 
                    tau_sq, 
                    omega_sq, 
                    bal, 
                    rho, 
                    P, 
                    f_c_val,
                    sigma_j_sq_inc = NULL,
                    pilot_data = NULL,
                    return_study_params = FALSE,
                    seed = NULL,
                    summarize_results = FALSE){
  
  #REMOVE LATER
  #start_time = Sys.time()
  
  # in this case we are specifying power unconditionally. If we were able to we would specify power unconditionally over the distribution of study features that are generated. 
  #but that is tricky to do since unconditional power will then depend on the whole distribution of study features. 
  # So, we need to simplify it a little for the purpose of the sim. still specify range of power of levels for a data set that has all average study features
  # meaning now variation in the number of effects or primary study sample sizes. 
  
  if (!is.null(seed)) set.seed(seed)
  
  
  # for determining the mu vector. but for k_j and sigma_j_sq and N use the average of empirical distribution. will result in fixed mu values. 
  # if we generate data that is pretty balanced then true power level should correspond exactly to what we specify. 
  # if kj and sigma are not constant in data generated, it will not map on exactly but should be in the neighborhood range power of levels. 
  
  # used average of the empirical distribution.
  # will hold the mu values fixed
  # balanced true level will correspond 
  # true power won't be exactly of what we specified. 
  # still have a range of power levels. 
  
  N_mean = mean(pilot_data$N) 
  k_mean = mean(pilot_data$kj)
  
  if(sigma_j_sq_inc == TRUE){
    sigma_j_sq_mean = mean(pilot_data$sigma_j_sq)
  } else{
    sigma_j_sq_mean = NULL
  }
  
  
   mu_vector <- mu_values(J = J, tau_sq = tau_sq, omega_sq = omega_sq, 
                          rho = rho, P = P,
                          f_c_val = f_c_val,
                          bal =bal, 
                          k_j = k_mean,  
                          sigma_j_sq = sigma_j_sq_mean, 
                          N = N_mean
                          )
  
  
   
   
  results <- map(1:iterations, ~{
    sample_dat <- n_ES_empirical(pilot_data, J = J)
    
 
    
    if(sigma_j_sq_inc == TRUE){
      sigma_j_sq = sample_dat$sigma_j_sq
    } else{
      sigma_j_sq = NULL
    }
   
    
    dat <- generate_meta(J = J, 
                         tau_sq = tau_sq, 
                         omega_sq = omega_sq, 
                         bal = bal, 
                         rho = rho, 
                         P = P, 
                         sample_sizes = sample_dat$N, 
                         k_j = sample_dat$kj,
                         mu_vector = mu_vector,
                     
                         f_c_val = f_c_val,
                         sigma_j_sq = sigma_j_sq,
                         return_study_params = return_study_params,
                         seed = seed ## maybe should remove this here and only keep in run_sim
    )
    
    
    est_res <-  estimate_model(data = dat,
                               moderator_val = dat$category,
                               cluster_id = dat$studyid,
                               delta = dat$g, 
                               delta_var =  dat$var_g,
                               es_id = dat$esid,
                               r= rho,
                               smooth_vi = TRUE, 
                             #  return_mu = NULL,
                               control_list = list()
    )
    
   
    
    est_res$mu_vector_list <-  list(mu_vector)
    
    est_res
    
  }) |> dplyr::bind_rows()
  
  #REMOVE LATER
  #end_time = Sys.time()
  
  #results$time <-   end_time - start_time
  
  if (summarize_results) {
    performance <- sim_performance(results = results)
    return(performance)
  } else {
    return(results)
  }
}

