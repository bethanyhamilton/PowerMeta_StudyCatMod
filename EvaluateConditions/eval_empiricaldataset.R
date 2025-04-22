# Math data set -- williams et al.
# need to be cautious of the big studies. studies are CRT --- may have small SEs.
# will find effective sample size rather than the raw sample size
# use the cluster adjusted vi and not N

#Note: 
#in first run: I actually averaged the vi to the 
#SampleID and summarized how many kj there were in each SampleID 
#(and not the CitationID). 
#This meta-analysis had 191 studies with 510 samples and 1109 effect 
#sizes. Then, after removing ESS >= 500 (where ESS = 4/(vi)) 
#and kj >= 20, the dataset becomes 178 studies and 463 samples. 
#Finally, after removing kj > N-2,  the dataset becomes 177 studies and 462 samples.

#While JEP thinks it's well-justified whether I average to the 
#study-level or to the sample-level to build the empirical datasets,
#I decided for the second run to average to the study-level instead to be more consistent
#with Vembye et al 2023. Doing this instead,
#after removing ESS >= 500 (where ESS = 4/(variance of the effect size)) and 
#kj >= 20, the dataset goes from 191 to 169 studies. 
#After removing kj > N-2,  the dataset becomes 167 studies. 


# -------------------------------------------------------------
# Empirical Data set -- cluster by study
# -------------------------------------------------------------

rm(list=ls())


math_dat <- read.csv("EvaluateConditions/data_191_RCTs.csv")

math_dat <- math_dat |>
  dplyr::select(
    CitationID,
    SampleID,
    EffectSizeID,
    yi,
    vi,
    L3_StuSamSiz,
    L3_StuSamSiz_T,
    L3_StuSamSiz_C,
    L3_TchSamSiz,
    L3_TchSamSiz_T,
    L3_TchSamSiz_C,
    L3_SchSamSiz,
    L3_SchSamSiz_T,
    L3_SchSamSiz_C
  )

math_dat <-
  math_dat |>
  mutate(ESS = 4 / vi#,
         # sample ID
       #  cluster_id = paste(CitationID, "_", SampleID, sep = "")
       ) |>
  group_by(CitationID) |>
  summarise(
 #   CitationID = unique(CitationID),
    g_avg = mean(yi),
    se_avg = mean(sqrt(vi)),
    kj = n(),
    N_ess = as.integer(round(mean(ESS))),
    # N_org = as.integer(round(mean(total_N))),
    # N_diff = N_org - N_ess,
    .groups = "drop"
  ) |>
  mutate(N_range = cut(N_ess, breaks = c(0, 500, 1000, 5000, 10000, 50000)))

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

ggplot(math_dat, aes(se_avg^2)) + 
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

math_dat |> 
  summarise(
 #   num_samples = length(unique(cluster_id)),
    num_studies = length(unique(CitationID))
  )


# only including studies where N < 500, and kj<20
dat_kjN_math <- 
  math_dat |>  
  dplyr::select(CitationID, kj, N = N_ess, se_avg) |> 
  mutate(sigma_j_sq = se_avg^2) |> 
  # Excluding effective sample sizes larger than 500 
  # and studies with more than 20 outcomes
  filter(N < 500, kj < 20) |> 
  as.data.frame()

dat_kjN_math |> 
  summarise(
    num_studies = length(unique(CitationID))
  )

# only including samples where kj <= N-2
dat_kjN_math <- dat_kjN_math |> filter(kj <= N - 2)

dat_kjN_math |> 
  summarise(
    num_studies = length(unique(CitationID))
  )

range(dat_kjN_math$kj)
range(dat_kjN_math$sigma_j_sq)

#Arithmetic means
mean(dat_kjN_math$kj)
mean(dat_kjN_math$sigma_j_sq)

# Harmonic means
dat_kjN_math |> 
  mutate(inv_kj = 1/kj,
         inv_sigma_j_sq = 1/sigma_j_sq) |> 
  summarise(kj_harmonic_mean = n()/sum(inv_kj),
            sigma_j_sq_harmonic_mean = n()/sum(inv_sigma_j_sq))



# Gamma parameters for N and sigma_j_sq
shape_rate1 <- MASS::fitdistr(dat_kjN_math$N, "gamma")
shape_rate2 <- MASS::fitdistr(dat_kjN_math$se_avg^2, "gamma")

library(fitdistrplus)
dist <- fitdistrplus::fitdist(dat_kjN_math$N, "gamma")
plot(dist)
dist <- fitdistrplus::fitdist(dat_kjN_math$sigma_j_sq, "gamma")
plot(dist)



# pmf for kj
set.seed(20250419)

lambda <- mean(dat_kjN_math$kj) - 1
x_val = 1:max(dat_kjN_math$kj)

stylized_dat <- data.frame(x = x_val, y = dpois(x_val - 1 ,  lambda), group = "Stylized")

empirical_dat <- as.data.frame(table(dat_kjN_math$kj)) |> mutate(x = as.integer(Var1),
                                                                 y = Freq/sum(Freq),
                                                                 group = "Empirical") |> dplyr::select(-Var1, -Freq)

all_dat <- bind_rows(stylized_dat, empirical_dat)


ggplot(all_dat, aes(x = x, y = y, fill = group)) +
  geom_col(position = "identity", alpha = 0.6, width = 0.4) +
  scale_fill_manual(values = c("Stylized" = "#F8766D", "Empirical" = "#00BFC4"),
                    labels = c(
                      "Empirical" = expression(paste("Empirical PMF for ", k[j])),
                      "Stylized" = expression(paste("Stylized PMF for ", k[j], " (Shifted Poisson)"))
                    )) +
  labs(x = expression(paste(k[j])), y = "Probability / Relative Frequency",
       title = expression(paste("PMF Plots of ", k [j], " from Epirical vs Stylized Sampling Methods")),
       fill = "") +
  theme_bw()



# density plot sigma_j_sq

set.seed(20250419)

shape_rate <- MASS::fitdistr(dat_kjN_math$se_avg^2, "gamma")

x <- seq(0.008, max(dat_kjN_math$sigma_j_sq), length.out = 462)

gamma_density <- dgamma(x, shape = shape_rate$estimate[1], 
                        rate = shape_rate$estimate[2])

density_sigma_sq <- density(dat_kjN_math$sigma_j_sq)
empirical_dat <- data.frame(x = density_sigma_sq$x, 
                            y = density_sigma_sq$y, group = "Empirical")

stylized_dat <- data.frame(x = x, y = gamma_density, group = "Stylized (Gamma)")
all_dat <- bind_rows(empirical_dat, stylized_dat)

ggplot(all_dat, aes(x = x, y = y, color = group)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = c("Empirical" = "#00BFC4", "Stylized (Gamma)" = "#F8766D"),
    labels = c(
      "Empirical" = expression(paste("Empirical Density for ", sigma[j]^2)),
      "Stylized (Gamma)" = expression(paste("Stylized Density for ", sigma[j]^2, " (Gamma)"))
    ),
    guide = guide_legend(override.aes = list(size = 2))
  ) +
  labs(
    x = expression(sigma[j]^2),
    y = "Density",
    color = "",
    title = expression(paste("Density plots of ", sigma[j]^2, " from Empirical vs Stylized Sampling Methods"))
  )  +
  theme_bw()



saveRDS(dat_kjN_math, "SimFunctions/dat_kjN_mathdat.rds")





## previous study:
# rm(list=ls())
# library(haven)
# library(tidyverse)
# 
# Diet_dat <-
#   read_dta("EvaluateConditions/RER_cluster.dta")
# 
# 
# 
# Diet_dat  <- Diet_dat  |>
#   mutate(
#     ESS = 4 / SE_g_adj^2
#   ) |>
#   filter(Follow_up == 0) |>
#   group_by(
#     Ref_nr, Authors
#   ) |>
#   summarise(
#     g_avg = mean(Effectsize_g_adj),
#     se_avg = mean(SE_g),
#     se_adj = mean(SE_g_adj),
#     kj = n(),
#     N_ess = as.integer(round(mean(ESS))),
#     N_org = as.integer(round(mean(N_total))),
#     N_diff = N_org - N_ess,
#     .groups = "drop"
#   ) |>
#   mutate(
#     N_range = cut(N_ess, breaks = c(0,500,1000,5000,10000,50000))
#   ); Diet_dat
# 
# Diet_dat |>
#   count(N_range)
# 
# ggplot(Diet_dat, aes(N_ess)) +
#   geom_histogram() +
#   scale_x_log10()
# 
# Diet_dat |>
#   filter(N_ess <= 500) |>
#   ggplot(aes(N_ess)) +
#   geom_histogram()
# 
# ggplot(Diet_dat, aes(kj)) +
#   geom_histogram()
# 
# ggplot(Diet_dat, aes(kj, N_ess)) +
#   geom_point() +
#   scale_y_log10()
# 
# # Creating dataset used for extracting empirical distributions of kj and N
# # for each study
# 
# dat_kjN <-
#   Diet_dat |>
#   select(kj, N = N_ess, se_adj) |>
#   mutate(sigma_j_sq = se_adj^2) |> 
#   # Excluding effective sample sizes larger than 500
#   # and studies with more than 20 outcomes
#   filter(N < 500, kj < 20) |>
#   as.data.frame()
# 
# range(Diet_dat$kj)
# length(unique(Diet_dat$Authors))
# 
# shape_rate <- MASS::fitdistr(dat_kjN$N, "gamma")
# shape_rate
# shape_rate <- MASS::fitdistr(dat_kjN$sigma_j_sq, "gamma")
# shape_rate
# saveRDS(dat_kjN, "EvaluateConditions/dat_kjN_Diet_dat.rds")
