#-------------------------------------------------------------------------------
#---------------------- MSQ for balanced c1-model ------------------------------
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

pi_c1_msq <- function(histdat,
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
        
        # Meansquares following Satterthwaite 1941
        ms_star_a <- ssq[1]/(I-1)
        ms_star_b <- ssq[1]/(J-1)
        ms_star_res <- ssq[1]/((I-1)*(J-1))
        
        #-----------------------------------------------------------------------
        
        # Gewichte
        w_a <- 1+1/I
        w_b <- 1+1/J
        w_res <- 1-1/(I*J)
        
        # df Satterthwaite
        
        df_s <- (sum(w_a*ms_star_a,
                    w_b*ms_star_b,
                    w_res*ms_star_res)^2) /
                sum((w_a*ms_star_a^2) / df["a"],
                    (w_b*ms_star_b^2) / df["b"],
                    (w_res*ms_star_res^2) / df["res"]
                    )
        
        # Please note,that the Statterthwaite approximation was done based on
        # the whole prediction variance var(y^*)=var(y)+var(mu) instead of 
        # the variance of the data var(y). Hence the approx. df 3.352 used for 
        # interval calculation differs from that used by Hoffman & Berger (6.96)
        
        #-----------------------------------------------------------------------
        
        # Gesamtvarianz (var_y + var_y_mean)
        
        var_total <- sum(w_a*ms_star_a,
                         w_b*ms_star_b,
                         w_res*ms_star_res)
                
        #-----------------------------------------------------------------------
        
        # PI for 1 future y
        y_mean <- mean(histdat[,3])
        L <- y_mean - qt(1-alpha/2, df=df_s)*sqrt(var_total)
        U <- y_mean + qt(1-alpha/2, df=df_s)*sqrt(var_total)
        pi_c2 <- c("L"=L, "U"=U)
        
        print(qt(1-alpha/2, df=df_s)*sqrt(var_total))
        #-----------------------------------------------------------------------
        
        output <- pi_c2
        
        return(output)
        
        
}


system.time(pi_msq <- pi_c1_msq(histdat=dat_hb))

# Interval on log scale
pi_msq

# interval on response scale
exp(pi_msq)



