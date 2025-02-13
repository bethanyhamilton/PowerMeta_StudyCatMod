rm(list=ls())
library(tidyverse)


source("SimFunctions/functionsdatagen.R")


# -------------------------------------------------------------
# Determining which combination of conditions are far-fetched
# -------------------------------------------------------------

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
    f_c_val = f_c_val, 
    sigma_j_sq = NULL,
    bal = "balanced_j"
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


#save(results, file = "EvaluateConditions/approximation_results.RData")

load("EvaluateConditions/approximation_results.RData")


View(results |> filter(f_c_val == 1))

View(results |> filter(f_c_val == 6))

results <- results|>
  rowwise() |> mutate(
  max_mu = max(as.numeric(unlist(mu_val))),
  J_C = J/C
) |> mutate(
  f_c_val_lab = case_when(
    f_c_val == 1 ~'f = [0, 1]',
    f_c_val == 2 ~'f = [0, 1, 2]',
    f_c_val == 3 ~'f = [0, 0, 1]',
    f_c_val == 4 ~'f = [0, 1, 1]',
    f_c_val == 5 ~'f = [0, 1, 2, 3]',
    f_c_val == 6 ~'f = [0, 0, 1, 2]',
    f_c_val == 7 ~'f = [0, 0, 0, 1]',
    f_c_val == 8 ~'f = [0, 0, 1, 1]',
  ),
  
  # f_c_val_lab = case_when(
  #   f_c_val == 1 ~'f = "["*0*", " *1*"]"',
  #   f_c_val == 2 ~'f = "["*0*", " *1*", " *2*"]"',
  #   f_c_val == 3 ~'f = "["*0*", " *0*", " *1*"]"',
  #   f_c_val == 4 ~'f = "["*0*", " *1*", " *1*"]"',
  #   f_c_val == 5 ~'f = "["*0*", " *1*", " *2*", " *3*"]"',
  #   f_c_val == 6 ~'f = "["*0*", " *0*", " *1*", " *2*"]"',
  #   f_c_val == 7 ~'f = "["*0*", " *0*", " *0*", " *1*"]"',
  #   f_c_val == 8 ~'f = "["*0*", " *0*", " *1*", " *1*"]"',
  # ),
  
  tau_sq = as.character(tau_sq),
  
  tau_lab = case_when(
    tau_sq == 0.0025 ~ "0.05",
    tau_sq == 0.04 ~ "0.20",
    tau_sq == 0.16 ~ "0.40",
    
  ),
  omega_sq = as.character(omega_sq),
  omega_lab = case_when(
    omega_sq == 0.0025 ~ "0.05",
    omega_sq == 0.01 ~ "0.10",
    omega_sq == 0.04 ~ "0.20",
    
  ),
  
  rho_lab = case_when(
    rho == 0.2 ~ ".2",
    rho == 0.5 ~ ".5",
    rho == 0.8 ~ ".8",
    
  ),
  
  J_lab = case_when(
    J == 12 ~ "J = 12",
    J == 24 ~ "J = 24",
    J == 36 ~ "J = 36",
    J == 48 ~ "J = 48",
    J == 60 ~ "J = 60",
    J == 72 ~ "J = 72",
    
  ),
  
  J_C_lab = case_when(
    J_C == 6 ~ "J_C = 6",
    J_C == 12 ~ "J_C = 12",
    J_C == 18 ~ "J_C = 18",
    J_C == 24 ~ "J_C = 24",
    J_C == 30 ~ "J_C = 30",
    J_C == 36 ~ "J_C = 36",
    J_C == 4 ~ "J_C = 4",
    J_C == 8 ~ "J_C = 8",
    J_C == 16 ~ "J_C = 16",
    J_C == 20 ~ "J_C = 20",
    J_C == 3 ~ "J_C = 3",
    J_C == 9 ~ "J_C = 9",
    J_C == 15 ~ "J_C = 15",
    
    
  )
  
 
  
  
  
)

results$J_C_lab <- factor(results$J_C_lab, levels = c("J_C = 3", "J_C = 4", "J_C = 6",
                                     "J_C = 8", "J_C = 9", "J_C = 12",
                                     "J_C = 15", "J_C = 16", "J_C = 18",
                                     "J_C = 20", "J_C = 24", "J_C = 30", 
                                     "J_C = 36"))

results$f_c_val_lab <- factor(results$f_c_val_lab, levels = c("f = [0, 1]" ,
                                                              "f = [0, 1, 2]" ,
                                                              "f = [0, 0, 1]",
                                                              "f = [0, 1, 1]"  ,
                                                              "f = [0, 1, 2, 3]",
                                                              "f = [0, 0, 1, 2]",
                                                              "f = [0, 0, 0, 1]",
                                                              "f = [0, 0, 1, 1]"
  
))


results$J_lab <- factor(results$J_lab, levels = c("J = 12", "J = 24",
                                                  "J = 36", "J = 48", 
                                                  "J = 60", "J = 72"))

max(results$max_mu)

### filter max value (1.5 SD) ** J= 12 condition probably. 
### reducing rho and omega values to two condition possibly.

#### Graphs

# pattern_labs <- c("f = [0, 1]", "f = [0, 1, 2]", "f = [0, 0, 1]", "f = [0, 1, 1]", "f = [0, 1, 2, 3]",
#                   "f = [0, 0, 1, 2]", "f = [0, 0, 0, 1]", "f = [0, 0, 1, 1]")
# #names(pattern_labs) <- c("1", "2", "3", "4", "5", "6", "7", "8")



#tau_sq_labs <- c('tau^2*" = 0.0025"', 'tau^2*" = 0.04"', 'tau^2*" = 0.16"')
#names(pattern_labs) <- c("1", "2", "3", "4", "5", "6", "7", "8")




### x = power, y= max mu, tau ~ pattern, color = j
results |>
  ggplot(aes(x = as.character(P), y = max_mu, 
             color = J_lab, shape = rho_lab)) +
  geom_point(alpha = .5, position = position_jitter(width = .3)) + 
  facet_grid(tau_lab~f_c_val_lab, labeller = label_bquote(rows = tau == . (tau_lab),
                                                          cols = .(f_c_val_lab))) +
  labs(
    x = "Power Level", 
    y = expression("Max " *mu* " Value"),
    color = "", shape = ""
  ) + 
  theme_bw() +
  theme(
    plot.caption=element_text(hjust = 0, size = 10),
    legend.position= "bottom",
    panel.spacing.x = unit(5, "mm"), 
    panel.spacing.y = unit(5, "mm"),
  )


results_test <- results |> filter( rho_lab == ".5") |> filter(J_lab  == "J = 36")


results_test |>
  ggplot(aes(x = as.character(P), y = max_mu, 
             color = omega_lab)) +
  geom_point(alpha = .5, position = position_jitter(width = .3)) + 
  facet_grid(tau_lab~f_c_val_lab, labeller = label_bquote(rows = tau == . (tau_lab),
                                                          cols = .(f_c_val_lab))) +
  labs(
    x = "Power Level", 
    y = expression("Max " *mu* " Value"),
    color = "", shape = ""
  ) + 
  theme_bw() +
  theme(
    plot.caption=element_text(hjust = 0, size = 10),
    legend.position= "bottom",
    panel.spacing.x = unit(5, "mm"), 
    panel.spacing.y = unit(5, "mm"),
  )

### x = power, y= max mu, tau ~ pattern, color = j_c
results |>
  ggplot(aes(x = as.character(P), y = max_mu,
             color = J_C_lab)) +
  geom_point(alpha = .5) +
  facet_grid(tau_lab~f_c_val_lab, labeller = label_bquote(rows = tau == . (tau_lab),
                                                          cols = .(f_c_val_lab))) +
  labs(
    x = "Power Level",
    y = expression("Max " *mu* " Value"),
    color = "", shape = ""
  ) +
  theme_bw() +
  theme(
    plot.caption=element_text(hjust = 0, size = 10),
    legend.position= "bottom",
    panel.spacing.x = unit(5, "mm"),
    panel.spacing.y = unit(5, "mm"),
  )



### x = power, y= max mu, rho ~ pattern, color = J
  results |> 
  ggplot(aes(x = as.character(P), y = max_mu, 
             color = J_lab, 
          
            )
         ) +

  geom_point(alpha = .5) + 
  facet_grid(rho_lab~f_c_val_lab, labeller = label_bquote(rows = rho == . (rho_lab),
                                                            cols = .(f_c_val_lab))) +

  labs(
    x = "Power Level", 
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

  ### x = power, y= max mu, tau ~ omega, color = J
  results |> 
    ggplot(aes(x = as.character(P), y = max_mu, 
               color = J_lab, 
           
    )
    ) +

    geom_point(alpha = .5) + 
    facet_grid(tau_lab~omega_lab, 
               labeller =label_bquote(rows = tau == .(tau_lab),
                                      cols = omega == .(omega_lab)
               ) ) +
    labs(
      x = "Power Level", 

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
  
  