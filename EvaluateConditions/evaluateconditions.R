rm(list=ls())
library(tidyverse)


source("SimFunctions/functionsdatagen.R")




run_mu <- function(
                   J, 
                   tau_sq, 
                   omega_sq, 
                   rho, 
                   P, 
                   k_j = 3, ####CHANGE THIS LATER
                   N = 30, ####CHANGE THIS LATER
                   f_c_val ) {
  
  mu_vector <- mu_values(
    J = J,
    tau_sq = tau_sq,
    omega_sq = omega_sq,
    rho = rho,
    P = P,
    k_j = k_j,
    N = N,
    f_c_val = f_c_val
  )
  
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
  
  
  return(
    tibble(
      C= C, 
      mu_val = list(mu_vector)
    )
  )
  
  
}





design_factors_bal <- list(
 # C = c(2, 3, 4),
  J = c(12, 24, 36, 48, 60, 72),
  tau_sq = c(0.05, 0.20, 0.40)^2, 
  omega_sq = c(0.05, 0.10, 0.20)^2,
  rho = c(.2, .5, .8),
  P = c(0.05, 0.2, 0.4, 0.6, 0.8, 0.9),
  f_c_val = c(1, 2, 3, 4, 5, 6, 7, 8)
)


params <- expand.grid(design_factors_bal) 
dim(params) # 7776, so about 3.888 min to run  

tm <- system.time(results <- plyr::mdply(params, 
.fun = run_mu,
.inform = TRUE,
.parallel = FALSE))


save(results, file = "EvaluateConditions/approximation_results.RData")


View(results |> filter(f_c_val == 1))

View(results |> filter(f_c_val == 6))

results <- results|>
  rowwise() |> mutate(
  max_mu = max(as.numeric(unlist(mu_val))),
  J_C = J/C
)


results %>% 
  ggplot(aes(x = as.character(P), y = max_mu, 
             color = as.character(J_C))) +
  geom_point(alpha = .5) + 
  facet_grid(tau_sq~f_c_val) +
  labs(
    x = "Approximated power", 
    y = expression("Max  " *mu* "  Value"),
    color = "", shape = ""
  ) + 
  theme_bw() +
  theme(
    plot.caption=element_text(hjust = 0, size = 10),
    legend.position= "bottom",
    panel.spacing.x = unit(5, "mm"), 
    panel.spacing.y = unit(5, "mm"),
  )

  results %>% 
  ggplot(aes(x = as.character(P), y = max_mu, 
             color = as.character(J_C), 
          
            )
         ) +

  geom_point(alpha = .5) + 
  facet_grid(rho~f_c_val) +
  labs(
    x = "Approximated power", 
   y = expression("Max  " *mu* "  Value"),
    color = "", shape = ""
  ) + 
  theme_bw() +
  theme(
    plot.caption=element_text(hjust = 0, size = 10),
    legend.position= "bottom",
    panel.spacing.x = unit(5, "mm"), 
    panel.spacing.y = unit(5, "mm"),
  )

  results %>% 
    ggplot(aes(x = as.character(P), y = max_mu, 
               color = as.character(J_C), 
           
    )
    ) +

    geom_point(alpha = .5) + 
    facet_grid(tau_sq~omega_sq) +
    labs(
      x = "Approximated power", 

      y = expression("Max  " *mu* "  Value"),
      color = "", shape = ""
    ) + 
    theme_bw() +
    theme(
      plot.caption=element_text(hjust = 0, size = 10),
      legend.position= "bottom",
      panel.spacing.x = unit(5, "mm"), 
      panel.spacing.y = unit(5, "mm"),
    )
  

