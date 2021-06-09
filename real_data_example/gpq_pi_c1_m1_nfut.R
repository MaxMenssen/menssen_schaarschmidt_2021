#-------------------------------------------------------------------------------
#---------------------- GPQ Lin & Liao Method 3 --------------------------------
#-------------------------------------------------------------------------------

# library(mvtnorm)

#-------------------------------------------------------------------------------
# Data of Table 4 from Hoffmann & Berger 2011 (Normalized mean RU)

run_1 <- c(1.11, 1.00, 1.50, 1.04, 0.929, 1.02, 0.965, 0.973, 1.01, 1.02,
           1.04, 1.13, 1.03, 1.04, 1.02, 1.04, 0.929, 1.04, 0.947, 0.947)

run_2 <- c(1.25, 1.18, 1.52, 1.15, 1.13, 1.03, 1.18, 1.18, 1.21, 1.29,
           1.16, 1.22, 1.35, 1.28, 1.29, 1.38, 1.20, 1.16, 1.20, 1.10)

run_3 <- c(1.03, 0.961, 1.30, 1.12, 1.06, 1.10, 0.804, 0.850, 1.02, 1.07, 
           0.933, 0.979, 1.11, 1.12, 1.02, 1.03, 0.915, 1.02, 0.915, 0.961)

run <- factor(c(rep(1, 20), rep(2, 20), rep(3, 20)))

mouse_id <- factor(c(1:20, 1:20, 1:20))

dat_hb <- data.frame(run=run,
                     mouse_id=mouse_id,
                     log_nmru=log(c(run_1, run_2, run_3)))

#-------------------------------------------------------------------------------

pi_c1_gpq_m3 <- function(histdat,
                         GPQ_n=1e+04,
                         alpha=0.05){
  
  
  #-----------------------------------------------------------------------------
  # number of model terms
  
  # I
  I <- length(levels(droplevels(histdat[,1])))
  
  # J
  J <- length(levels(droplevels(histdat[,2])))
  
  #-----------------------------------------------------------------------------
  
  # modelfit
  c1_aov <- aov(histdat[,3] ~ histdat[,1] + histdat[,2], data=histdat)
  
  #-----------------------------------------------------------------------------
  
  # degrees of freedom
  df <- summary(c1_aov)[[1]][["Df"]]
  names(df) <- c("a", "b", "res")
  
  #-----------------------------------------------------------------------------
  # sum of sqares
  ssq <- summary(c1_aov)[[1]][["Sum Sq"]]
  names(ssq) <- c("a", "b", "res")
  
  #-------------------------------------------------------------------------------
  
  # GPQ´s for Psi_Squares (step 3)
  R_Psi_a <- ssq["a"] / rchisq(n=GPQ_n, df=df["a"])
  R_Psi_b <- ssq["b"] / rchisq(n=GPQ_n, df=df["b"])
  R_Psi_e <- ssq["res"] / rchisq(n=GPQ_n, df=df["res"])
  
  #-----------------------------------------------------------------------------
  
  # GPQ for variance components (step 4)
  
  # ATTENTION: Some estimates for the variance components can become negative
  # and are adjusted to be at least 0
  
  R_sig2_a_v <- 1/(J) * (R_Psi_a-R_Psi_b)
  R_sig2_b_v <- 1/(I) * (R_Psi_b-R_Psi_b)
  R_sig2_e_v <- R_Psi_e
  
  R_sig2_a <- as.list(replace(R_sig2_a_v, R_sig2_a_v<0, 0))
  R_sig2_b <- as.list(replace(R_sig2_b_v, R_sig2_b_v<0, 0))
  R_sig2_e <- as.list(replace(R_sig2_e_v, R_sig2_e_v<0, 0))
  
  R_sig2_y <- (Map("+", Map("+", R_sig2_a, R_sig2_b), R_sig2_e))
  
  #-----------------------------------------------------------------------------
  
  # GPQ for var(mu) (step 5b)
  GPQ_sig2_mu <- as.list(1/(I*J) * (J*as.numeric(R_sig2_a) +
                                      I*as.numeric(R_sig2_b) + 
                                      as.numeric(R_sig2_e)))
  
  #-----------------------------------------------------------------------------
  
  # GPQ for var(y_star) (step 5)
  GPQ_sigma_k <- Map("+", R_sig2_y, GPQ_sig2_mu)
  
  # Esimating the prediction error
  D_est_k <- qnorm(p=1-alpha/2, mean=0, sd=sqrt(as.numeric(GPQ_sigma_k)))
  
  D_est <- median(D_est_k)
  
  #-----------------------------------------------------------------------------
  
  # PI for 1 future y
  y_mean <- mean(histdat[,3])
  L <- y_mean - D_est
  U <- y_mean + D_est
  pi_c2 <- c("L"=L, "U"=U, "D_est"=D_est)
  
  #----------------------------------------------------------------------------
  
  # Output
  output <- pi_c2
  
  return(output)
}

#------------------------------------------------------------------------------

# Interval calculation
system.time(pi_gpq_m3 <- pi_c1_gpq_m3(histdat=dat_hb))

# Interval on log scale
pi_gpq_m3

# interval on response scale
exp(pi_gpq_m3)
