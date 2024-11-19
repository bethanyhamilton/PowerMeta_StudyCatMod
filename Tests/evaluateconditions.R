source("../SimFunctions/functionsdatagen.R")


design_factors_bal <- list(
  C = c(2, 3, 4),
  #C = c(2, 3, 4, 5, 6),
  #f_c = c(),
  J_c = c(3, 5, 10, 20),
  tau = c(0.05, 0.20, 0.40), 
  omega = c(0.05, 0.10, 0.20),
  rho = c(.2, .5, .8),
  P = c(0.05, 0.2, 0.4, 0.6, 0.8, .9)
)


## 3 categories, J_c = 3
data_ex <- tibble(studyid = c(1:9),
                  k_j = rep(3, 9),
                  n_j = rep(30, 9),
                  sigma_j = sqrt(4 / n_j), 
                  omega = rep(.10, 9),
                  rho = rep(.5, 9),
                  tau = rep(.20, 9),
                  three_cat =c(rep("a",3), rep("b",3), rep("c",3))
)


df2 <- multiple_categories(dat = data_ex, moderator = data_ex$three_cat, cluster = data_ex$studyid, sigma_j = data_ex$sigma_j , rho = data_ex$rho, omega_sq = data_ex$omega, tau_sq = data_ex$tau)

df_num <- df2$df_num[1]
df_den <- df2$nu_D[1] - df2$df_num[1] + 1

f_crit <- qf(1-.05, df_num,df_den )

lambda <- find_lambda(d1= df_num, d2 = df_den, x = f_crit, area = 1 - .9, interval = c(0, 100), tol = 0.0001)

f_c <- c(0, 1, 2)
#f_c <- c(0, 0, 1)
#f_c <- c(0, 1, 1)

zeta_val <- zeta(pattern = f_c, lambda = lambda, weights =df2$W )

build_mu(pattern = f_c, zeta = zeta_val)


