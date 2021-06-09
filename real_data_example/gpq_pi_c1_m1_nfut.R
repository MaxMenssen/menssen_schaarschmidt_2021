#-------------------------------------------------------------------------------
#---------------------- GPQ Lin & Liao Method 1 --------------------------------
#-------------------------------------------------------------------------------


source("c2.R" )
library(mvtnorm)


#-------------------------------------------------------------------------------

pi_c2_gpq <- function(mu, n_i, n_j, n_ij,
                      var_i, var_j, var_ij, var_ijk,
                      n_i_fut, n_j_fut, n_ij_fut,
                      # GPQ_n=1e+04,
                      GPQ_n=100,
                      alpha=0.05){
  
  
  sd_i <- sqrt(var_i)
  sd_j <- sqrt(var_j)
  sd_ij <- sqrt(var_ij)
  sd_ijk <- sqrt(var_ijk)
  
  # Historical data from h1 design
  c2dat <-interblock(mu=mu, n_i=n_i, n_j=n_j, n_ij=n_ij,
                     sd_i=sd_i, sd_j=sd_j, sd_ij=sd_ij, sd_ijk=sd_ijk)
  
  #-----------------------------------------------------------------------------
  # number of model terms
  
  # I
  I <- length(levels(droplevels(c2dat$a)))
  
  # J
  J <- length(levels(droplevels(c2dat$b)))
  
  # K
  K <- length(c2dat$a : c2dat$b) / length(levels(droplevels(c2dat$a : c2dat$b)))
  
  # IJ
  IJ <- length(levels(droplevels(c2dat$a : c2dat$b)))
  
  # IJK
  IJK <- nrow(c2dat)
  
  #-----------------------------------------------------------------------------
  # modelfit
  
  c2_aov <- aov(y_ijk ~ a*b, c2dat)
  
  #-----------------------------------------------------------------------------
  
  # degrees of freedom
  df <- summary(c2_aov)[[1]][["Df"]]
  names(df) <- c(attributes(c2_aov$terms)$term.labels, "res")
  df
  
  #-----------------------------------------------------------------------------
  # sum of sqares
  ssq <- summary(c2_aov)[[1]][["Sum Sq"]]
  names(ssq) <- c(attributes(c2_aov$terms)$term.labels, "res")
  # ssq
  
  #-----------------------------------------------------------------------------
  # mean squares
  msq <- summary(c2_aov)[[1]][["Mean Sq"]]
  names(msq) <- c(attributes(c2_aov$terms)$term.labels, "res")
  
  
  #-----------------------------------------------------------------------------
  
  I_i <- diag(n_i_fut)
  J_i <- matrix(1,n_i_fut,n_i_fut)
  
  I_j <- diag(n_j_fut)
  J_j <- matrix(1,n_j_fut,n_j_fut)
  
  I_k <- diag(n_ij_fut)
  J_k <- matrix(1,n_ij_fut,n_ij_fut)
  
  #-------------------------------------------------------------------------------
  
  # GPQ´s for Psi_Squares (step 3)
  R_Psi_a <- ssq["a"] / rchisq(n=GPQ_n, df=df["a"])
  R_Psi_b <- ssq["b"] / rchisq(n=GPQ_n, df=df["b"])
  R_Psi_ab <- ssq["a:b"] / rchisq(n=GPQ_n, df=df["a:b"])
  R_Psi_e <- ssq["res"] / rchisq(n=GPQ_n, df=df["res"])
  
  #-----------------------------------------------------------------------------
  
  # GPQ for variance components (step 4)
  
  # ATTENTION: Some estimates for the variance components can become negative
  # and are adjusted to be at least 0
  
  R_sig2_a_v <- 1/(J*K) * (R_Psi_a-R_Psi_ab)
  R_sig2_b_v <- 1/(I*K) * (R_Psi_b-R_Psi_ab)
  R_sig2_ab_v <- 1/(K) * (R_Psi_ab-R_Psi_e)
  R_sig2_e_v <- R_Psi_e
  
  R_sig2_a <- as.list(replace(R_sig2_a_v, R_sig2_a_v<0, 0))
  R_sig2_b <- as.list(replace(R_sig2_b_v, R_sig2_b_v<0, 0))
  R_sig2_ab <- as.list(replace(R_sig2_ab_v, R_sig2_ab_v<0, 0))
  R_sig2_e <- as.list(replace(R_sig2_e_v, R_sig2_e_v<0, 0))
  
  #-----------------------------------------------------------------------------
  
  # GPQ for var(y) (step 5a)
  f_a <- function(X){X*(I_i %x% J_j %x% J_k)}
  f_b <- function(X){X*(J_i %x% I_j %x% J_k)}
  f_ab <- function(X){X*(I_i %x% I_j %x% J_k)}
  f_e <- function(X){X*(I_i %x% I_j %x% I_k)}
  
  vc_a <- lapply(X=R_sig2_a, FUN=f_a)
  vc_b <- lapply(X=R_sig2_b, FUN=f_b)
  vc_ab <- lapply(X=R_sig2_ab, FUN=f_ab)
  vc_e <- lapply(X=R_sig2_e, FUN=f_e)
  
  sigma_star_1 <- Map("+",vc_a, vc_b)
  sigma_star_2 <- Map("+", sigma_star_1, vc_ab)
  sigma_star <- Map("+", sigma_star_2, vc_e)
  
  #-----------------------------------------------------------------------------
  
  # GPQ for var(mu) (step 5b)
  GPQ_sig2_mu <- 1/IJK * (J*K*as.numeric(R_sig2_a) +
                            I*K*as.numeric(R_sig2_b) + 
                            K*as.numeric(R_sig2_ab) +
                            as.numeric(R_sig2_e))
  
  J_mu <- matrix(1, (n_i_fut*n_j_fut*n_ij_fut), (n_i_fut*n_j_fut*n_ij_fut))
  
  f_J_mu <- function(X){X*J_mu}
  
  sig2_mu_J_mu <- lapply(GPQ_sig2_mu, f_J_mu)
  
  # print(sig2_mu_J_mu)
  
  #-----------------------------------------------------------------------------
  
  # GPQ for var(y_star) (step 5)
  GPQ_sigma <- Map("+", sigma_star, sig2_mu_J_mu)
  
  # print(GPQ_sigma)
  
  # GPQ_var (step 6)
  mvn_quant <- function(X){qmvnorm(p=1-alpha/2, mean=0, sigma=X)$quantile}
  
  GPQ_D <- lapply(GPQ_sigma, mvn_quant)
  
  # print(as.numeric(GPQ_D))
  
  # GPQ_D <- qmvnorm(p=1-alpha/2, mean=0, sigma=GPQ_sigma)
  D_est <- median(as.numeric(GPQ_D))
  
  # print(D_est)
  
  #-----------------------------------------------------------------------
  
  
  # PI for 1 future y
  y_mean <- mean(c2dat$y_ijk)
  L <- y_mean - D_est
  U <- y_mean + D_est
  pi_c2 <- c("L"=L, "U"=U)
  
  
  #-----------------------------------------------------------------------
  
  # Future data
  c2dat_fut <- interblock(mu=mu, n_i=n_i_fut, n_j=n_j_fut, n_ij=n_ij_fut,
                          sd_i=sd_i, sd_j=sd_j, sd_ij=sd_ij, sd_ijk=sd_ijk)
  
  y_fut <- c2dat_fut$y_ijk
  
  # Coverage
  fut_in_pi_est <- all(all(L<y_fut), all(y_fut<U))
  
  
  # Output
  output <- list("pi_c2"=pi_c2,
                 "y_fut"=c2dat_fut,
                 "fut_in_pi_est"=fut_in_pi_est)
  
  return(output)
  
  
}

system.time(pi_est <- pi_c2_gpq(mu=80, n_i=4, n_j=3, n_ij=2,
                                n_i_fut=2, n_j_fut=2, n_ij_fut=2,
                                var_i=3, var_j=3, var_ij=3, var_ijk=3))
pi_est

# mit GPQ_n=100 dauert es 12.10 sec

# mit GPQ_n=10000 würde es es 20.16667 Minuten dauern
12.10*100/60


#-------------------------------------------------------------------------------


c2_coverage <- function(nsim, mu, n_i, n_j, n_ij,
                        n_i_fut, n_j_fut, n_ij_fut,
                        var_i, var_j, var_ij, var_ijk){
  
  futinint <- vector(length=nsim)
  
  for(i in 1:nsim){
    futinint[i] <- unlist(pi_c2_gpq(mu=mu, n_i=n_i, n_j=n_j, n_ij=n_ij,
                                    n_i_fut=n_i_fut, n_j_fut=n_j_fut, n_ij_fut=n_ij_fut,
                                    var_i=var_i, var_j=var_j, var_ij=var_ij, 
                                    var_ijk=var_ijk)["fut_in_pi_est"])
  }
  
  
  cover_prob <- sum(futinint)/nsim
  
  # defining the output object
  out <- c(nsim=nsim,
           mu=mu, 
           n_i=n_i, n_j=n_j, n_ij=n_ij,
           var_i=var_i, var_j=var_j, var_ij=var_ij, var_ijk=var_ijk, 
           cover_prob=cover_prob
  )
  
  return(out)
}



#-------------------------------------------------------------------------------


# Simulationsettings 
# simdat <- expand.grid(mu=c(80),
#                       I=c(5, 10, 15),
#                       J=c(2, 5, 10),
#                       IJ=c(2, 10),
#                       I_fut=2,
#                       J_fut=2, 
#                       IJ_fut=2,
#                       varI=2*c(10, 1,  0.1),
#                       varJ=2*c(10, 1, 0.1),
#                       varIJ=2*c(10, 1, 0.1),
#                       varIJK=2)

simdat <- expand.grid(mu=c(80),
                      I=c(5,  15),
                      J=c(2,  10),
                      IJ=c(2, 10),
                      I_fut=2,
                      J_fut=2, 
                      IJ_fut=2,
                      varI=2*c(10,   0.1),
                      varJ=2*c(10,  0.1),
                      varIJ=2*c(10,  0.1),
                      varIJK=2)

str(simdat)
#--------------------------------------------------------------------------------


system.time(
  sim_c2 <- apply(simdat, MARGIN=1,
                  FUN=function(X){
                    c2_coverage(nsim=1000,
                                mu=X[["mu"]],
                                n_i=X[["I"]],
                                n_j=X[["J"]],
                                n_ij=X[["IJ"]],
                                n_i_fut=X[["I_fut"]],
                                n_j_fut=X[["J_fut"]],
                                n_ij_fut=X[["IJ_fut"]],
                                var_i=X[["varI"]],
                                var_j=X[["varJ"]],
                                var_ij=X[["varIJ"]],
                                var_ijk=X[["varIJK"]])
                  })
)

#--------------------------------------------------------------------------------


system.time(
  sim_c2 <- apply(simdat, MARGIN=1,
                  FUN=function(X){
                    c2_coverage(nsim=5000,
                                mu=X[["mu"]],
                                n_i=X[["I"]],
                                n_j=X[["J"]],
                                n_ij=X[["IJ"]],
                                var_i=X[["varI"]],
                                var_j=X[["varJ"]],
                                var_ij=X[["varIJ"]],
                                var_ijk=X[["varIJK"]])
                  })
)

# simulation time
1616.90 *10 /60/60

results_c2 <- data.frame(t(sim_c2))

results_c2$method <- "gpq_m1"

write.csv2(results_c2, "gpq_c2_m1.csv")


