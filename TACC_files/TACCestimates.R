library(tidyverse)
library(metafor)
library(clubSandwich)
library(mvtnorm)

rm(list=ls())

source("SimFunctions/functionsdatagen.R")

design_factors <- list(
  J = c(24, 36, 48, 60, 72),
  tau_sq = c(0.05,  0.40)^2, 
  omega_sq = c(0.05,  0.20)^2,
  rho = c(.2,  .8),
  P = c(0.05, 0.2, 0.4, 0.6, 0.8, 0.9),
  f_c_val = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"),
  bal = c("balanced_j", "unbalanced_j") #,
#  batch_id = c("batch1","batch2") ### add this in run figure out new time for each batch when dividing 2500 replications across batches 
)

#run smaller batches of replications. set iterations to smaller so each row runs faster.
#instead of running in parallel we will run one condition/job on one node. This way when a job finishes in a queue another job will start
#this is more efficient than running in parallel in R.
#First pass do just batch id 1s. 
#maybe 1250 replications for first batch might be good -- **** 

# use pylaunch to run the everything on TACC

#save raw files



params <- expand.grid(design_factors) %>% 
  mutate(
    iterations = 20,
    seed = round(runif(1) * 2^30) + 1:n()
  ) %>%
  as_tibble()

set.seed(0303202571)
 params <- params |> 
   group_by(J) |> 
   sample_n(size=10, replace = TRUE) 



#  empirical_dat <- readRDS("SimFunctions/dat_kjN_mathdat.rds")
#  
#  
# 
# tm <- system.time(  results <- plyr::mdply(params,
#                          .fun = run_sim,
#                          pilot_data = empirical_dat,
#                          sigma_j_sq_inc = FALSE,
#                          .parallel = FALSE))
# 
# sum(results$time)
# tm
# 
# results <- list(tm, results)
# 
# saveRDS(results, tm, "SimFunctions/TACCtest.rds")

results <- readRDS("TACC_files/TACCtest.rds")

tm_j<- results |> 
  group_by(J) |> 
  summarize(mean_tm_20rep = mean(time),
            avg_tm_per_cond = mean_tm_20rep*(2500/20),
            tot_tm_byJ = avg_tm_per_cond*768) 

tm_j 

tm_all <- results |> ungroup() |> 
  summarise(mean_tm_20 = mean(time), 
            avg_tm_percond_all = mean_tm_20*(2500/20),
            avg_tm_percond_batch = mean_tm_20*(625/20),
            mean_tm_1rep = mean_tm_20/20)
tm_all


num_conditions <- 3840
num_core_pernode <- 48

## whole simulation


hrs_percond<- as.numeric((tm_all$avg_tm_percond_all))/(60^2)

hrs_acrosscond <-hrs_percond*num_conditions

hrs_per_node <- hrs_acrosscond/num_core_pernode

#across 4 nodes?
hrs_per_node/4


## with four batches of 625 
#625 rep per batch
hrs_perbatch<- as.numeric((tm_all$avg_tm_percond_batch))/(60^2)
hrs_perbatch_acrosscond <- hrs_perbatch*num_conditions
hrs_perbatch_node <- hrs_perbatch_acrosscond/num_core_pernode

#across 4 nodes?
hrs_perbatch_node/4



## actual time for 3840 conditions/jobs using the pylauncher 
## 4 batches of 3840 jobs each had 625 iterations 
## each batch was run on 4 nodes with 48 cores

batch1_min <- (3*60) + 30 + (47/60*1)
batch2_min <- (3*60) + 30 + (42/60*1)
batch3_min <- (3*60) + 31 + (47/60*1)
batch4_min <- (3*60) + 31 + (56/60*1)

one_job_min <- mean(batch1_min/3840, batch2_min/3840, 
                    batch3_min/3840, batch4_min/3840)


### if increase the number of conditions to 8640
## time in hrs for 8640 jobs four nodes
tm_8640jobs_hr_one_batch <- (one_job_min*8640)/60
#7.904375 hours for each batch on four nodes. 



