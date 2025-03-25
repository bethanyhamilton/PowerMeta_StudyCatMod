# -------------------------------------------------------------
# Empirical Data set 
# -------------------------------------------------------------

rm(list=ls())
library(haven)


Diet_dat <-
  read_dta("EvaluateConditions/RER_cluster.dta")



Diet_dat  <- Diet_dat  |>
  mutate(
    ESS = 4 / SE_g_adj^2
  ) |>
  filter(Follow_up == 0) |>
  group_by(
    Ref_nr, Authors
  ) |>
  summarise(
    g_avg = mean(Effectsize_g_adj),
    se_avg = mean(SE_g),
    se_adj = mean(SE_g_adj),
    kj = n(),
    N_ess = as.integer(round(mean(ESS))),
    N_org = as.integer(round(mean(N_total))),
    N_diff = N_org - N_ess,
    .groups = "drop"
  ) |>
  mutate(
    N_range = cut(N_ess, breaks = c(0,500,1000,5000,10000,50000))
  ); Diet_dat

Diet_dat |>
  count(N_range)

ggplot(Diet_dat, aes(N_ess)) +
  geom_histogram() +
  scale_x_log10()

Diet_dat |>
  filter(N_ess <= 500) |>
  ggplot(aes(N_ess)) +
  geom_histogram()

ggplot(Diet_dat, aes(kj)) +
  geom_histogram()

ggplot(Diet_dat, aes(kj, N_ess)) +
  geom_point() +
  scale_y_log10()

# Creating dataset used for extracting empirical distributions of kj and N
# for each study

dat_kjN <-
  Diet_dat |>
  select(kj, N = N_ess, se_adj) |>
  mutate(sigma_j_sq = se_adj^2) |> 
  # Excluding effective sample sizes larger than 500
  # and studies with more than 20 outcomes
  filter(N < 500, kj < 20) |>
  as.data.frame()

range(Diet_dat$kj)
length(unique(Diet_dat$Authors))

shape_rate <- MASS::fitdistr(dat_kjN$N, "gamma")
shape_rate
shape_rate <- MASS::fitdistr(dat_kjN$sigma_j_sq, "gamma")
shape_rate
saveRDS(dat_kjN, "EvaluateConditions/dat_kjN_Diet_dat.rds")
# #_______________________________________________________________________  
#   load("SimFunctions/refdata_igrm.Rdata")
#   
#   
#   ref_data <- 
#     datafull_igrm2 |> 
#     mutate(
#       ESS = 4 / vi_igrm
#     ) |> 
#     group_by(
#       studyID
#     ) |> 
#     summarise(
#       g_avg = mean(yi_igrm),
#       se_avg = mean(sqrt(vi_igrm)),
#       kj = n(),
#       N_ess = as.integer(round(mean(ESS))),
#       N_org = as.integer(round(mean(total_N))),
#       N_diff = N_org - N_ess,
#       .groups = "drop"
#     ) |> 
#     mutate(
#       N_range = cut(N_ess, breaks = c(0,500,1000,5000,10000,50000))
#     )
#   
#   ref_data |> 
#     count(N_range)
#   
#   ggplot(ref_data, aes(N_ess)) + 
#     geom_histogram() + 
#     scale_x_log10()
#   
#   ref_data |> 
#     filter(N_ess <= 500) |> 
#     ggplot(aes(N_ess)) + 
#     geom_histogram()
#   
#   ggplot(ref_data, aes(kj)) + 
#     geom_histogram()
#   
#   ggplot(ref_data, aes(kj, N_ess)) + 
#     geom_point() + 
#     scale_y_log10()
#   
#   # Creating dataset used for extracting empirical distributions of nj and N 
#   # for each study 
#   
#   dat_njN <- 
#     ref_data |>  
#     select(kj, N = N_ess) |> 
#     # Excluding effective sample sizes larger than 500 
#     # and studies with more than 20 outcomes
#     filter(N < 500, kj < 20) |> 
#     as.data.frame()
#   
#   length(unique(ref_data$studyID))
#   
#   range(ref_data$kj)
#________________________________________________________________________

erika_dat <- read.csv("EvaluateConditions/final_dataset.csv")



erika_dat <- erika_dat |> mutate(
  
  
  ES12..Control.Post.test.N = ifelse(is.na(ES12..Control.Post.test.N), ES19..Control.N, ES12..Control.Post.test.N),
  ES23..Intervention.N = ifelse(is.na(ES23..Intervention.N), ES20..Intervention.N, ES23..Intervention.N),
  ES12..Control.Post.test.N = ifelse(is.na(ES12..Control.Post.test.N), ES9..Total.sample.size/2, ES12..Control.Post.test.N),
  ES23..Intervention.N = ifelse(is.na(ES23..Intervention.N), ES9..Total.sample.size/2, ES23..Intervention.N),
  total_N = ES12..Control.Post.test.N + ES23..Intervention.N
  
)




erika_dat <- 
  erika_dat |> 
  mutate(
    ESS = 4 / igrm_d_v
  ) |> 
  group_by(
    StudyID
  ) |> 
  summarise(
    g_avg = mean(igrm_d),
    se_avg = mean(sqrt(igrm_d_v)),
    kj = n(),
    N_ess = as.integer(round(mean(ESS))),
    N_org = as.integer(round(mean(total_N))),
    N_diff = N_org - N_ess,
    .groups = "drop"
  ) |> 
  mutate(
    N_range = cut(N_ess, breaks = c(0,500,1000,5000,10000,50000))
  )

erika_dat |> 
  count(N_range)

ggplot(erika_dat, aes(N_ess)) + 
  geom_histogram() + 
  scale_x_log10()

erika_dat |> 
  filter(N_ess <= 500) |> 
  ggplot(aes(N_ess)) + 
  geom_histogram()

ggplot(erika_dat, aes(kj)) + 
  geom_histogram()

ggplot(erika_dat, aes(se_avg)) + 
  geom_histogram()

ggplot(erika_dat, aes(kj, N_ess)) + 
  geom_point() + 
  scale_y_log10()

ggplot(erika_dat, aes(kj, se_avg)) + 
  geom_point() + 
  scale_y_log10()


# Creating dataset used for extracting empirical distributions of nj and N 
# for each study 

dat_kjN <- 
  erika_dat |>  
  select(kj, N = N_ess, se_avg) |> 
  mutate(sigma_j_sq = se_avg^2) |> 
  # Excluding effective sample sizes larger than 500 
  # and studies with more than 20 outcomes
  filter(N < 500, kj < 20) |> 
  as.data.frame()

range(erika_dat$kj)

length(unique(dat_kjN$N))

saveRDS(dat_kjN, "EvaluateConditions/dat_kjN_erikadat.rds")

shape_rate <- MASS::fitdistr(dat_kjN$N, "gamma")
shape_rate
shape_rate <- MASS::fitdistr(dat_kjN$se_avg^2, "gamma")
shape_rate

# sum(is.na(erika_dat$total_N))


## NOTE make sure empirical dat has reports Var.

#__________________________________________________________________________

#Math data set

## no total N? I may not be able to use this one... 


# math -- williams. be cautious of the big studies. studies are CRT --- may have small SEs.
# need effective sample size rather than the raw sample size
# use the cluster adjusted sigma_j_sq ****and not N

rm(list=ls())


math_dat <- read.csv("EvaluateConditions/data_191_RCTs.csv")

math_dat <- math_dat |> select(CitationID, SampleID, EffectSizeID,yi, vi, L3_StuSamSiz, L3_StuSamSiz_T, L3_StuSamSiz_C, L3_TchSamSiz, L3_TchSamSiz_T, L3_TchSamSiz_C, L3_SchSamSiz, L3_SchSamSiz_T, L3_SchSamSiz_C )


math_dat <- 
  math_dat |> 
  mutate(
    ESS = 4 / vi,
    cluster_id = paste(CitationID, "_", SampleID, sep = "")
  ) |> 
  group_by(
    cluster_id
  ) |> 
  summarise(
    g_avg = mean(yi),
    se_avg = mean(sqrt(vi)),
    kj = n(),
    N_ess = as.integer(round(mean(ESS))),
   # N_org = as.integer(round(mean(total_N))),
   # N_diff = N_org - N_ess,
    .groups = "drop"
  ) |> 
  mutate(
    N_range = cut(N_ess, breaks = c(0,500,1000,5000,10000,50000))
  )

math_dat |> 
  count(N_range)

ggplot(math_dat, aes(N_ess)) + 
  geom_histogram() + 
  scale_x_log10()

math_dat |> 
  filter(N_ess <= 500) |> 
  ggplot(aes(N_ess)) + 
  geom_histogram()

ggplot(math_dat, aes(kj)) + 
  geom_histogram()

ggplot(math_dat, aes(se_avg)) + 
  geom_histogram()

ggplot(math_dat, aes(kj, N_ess)) + 
  geom_point() + 
  scale_y_log10()

ggplot(math_dat, aes(kj, se_avg)) + 
  geom_point() + 
  scale_y_log10()


## minimal SE -- something in neighborhood of this. 
2/sqrt(500) 

# Creating dataset used for extracting empirical distributions of nj and N 
# for each study 

dat_kjN_math <- 
  math_dat |>  
  select(kj, N = N_ess, se_avg) |> 
  mutate(sigma_j_sq = se_avg^2) |> 
  # Excluding effective sample sizes larger than 500 
  # and studies with more than 20 outcomes
  filter(N < 500, kj < 20) |> 
  as.data.frame()

dat_kjN_math<- dat_kjN_math |> filter(kj <= N - 2)

range(math_dat$kj)

length(unique(dat_kjN_math$N))

saveRDS(dat_kjN_math, "SimFunctions/dat_kjN_mathdat.rds")

shape_rate <- MASS::fitdistr(dat_kjN_math$N, "gamma")

shape_rate <- MASS::fitdistr(dat_kjN_math$se_avg^2, "gamma")

library(fitdistrplus)
dist <- fitdistrplus::fitdist(dat_kjN_math$N, "gamma")
plot(dist)

dist <- fitdistrplus::fitdist(dat_kjN_math$sigma_j_sq, "gamma")
plot(dist)


# density plot of se^2 and gamma distribution. 
# qqplot 
