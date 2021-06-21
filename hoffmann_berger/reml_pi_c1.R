#-------------------------------------------------------------------------------
#------ REML based PI with df approximation following van den Heuvel 2010-------
#-------------------------------------------------------------------------------

# Data of Table 4 from Hoffmann & Berger 2011 (Normalized mean RU)

# nmru for Run 1
run_1 <- c(1.11, 1.00, 1.50, 1.04, 0.929, 1.02, 0.965, 0.973, 1.01, 1.02,
           1.04, 1.13, 1.03, 1.04, 1.02, 1.04, 0.929, 1.04, 0.947, 0.947)

# nmru for Run 2
run_2 <- c(1.25, 1.18, 1.52, 1.15, 1.13, 1.03, 1.18, 1.18, 1.21, 1.29,
           1.16, 1.22, 1.35, 1.28, 1.29, 1.38, 1.20, 1.16, 1.20, 1.10)

# nmru for Run 3
run_3 <- c(1.03, 0.961, 1.30, 1.12, 1.06, 1.10, 0.804, 0.850, 1.02, 1.07, 
           0.933, 0.979, 1.11, 1.12, 1.02, 1.03, 0.915, 1.02, 0.915, 0.961)

# Variables for run and mouse_id
run <- factor(c(rep(1, 20), rep(2, 20), rep(3, 20)))
mouse_id <- factor(c(1:20, 1:20, 1:20))

# Attention: Ln(nmru)
dat_hb <- data.frame(run=run,
                     mouse_id=mouse_id,
                     log_nmru=log(c(run_1, run_2, run_3)))

#------------------------------------------------------------------------------

# PI following Franq et al. 

library(lme4)
library(VCA)

# Small helper function
nfun <- function(input){
        
        indvec <- input$re.assign$ind
        
        termsvec <- input$re.assign$terms
        
        termsind <- seq(from=1, to=length(termsvec))
        names(termsind) <- termsvec
        
        
        nvec <- numeric(length=length(termsind))
        names(nvec) <- termsvec
        
        for(i in 1:length(termsind)){
                
                nvec[i] <- sum(termsind[i]==indvec)
                
        }
        
        res <- nrow(input$data)
        names(res) <- "error"
        
        out <- c(nvec, res)
        
        return(out)
        
}

remlPI <- function(modelfit, type="both", alpha=0.05, printdfc=TRUE){
        
        # Variance components
        varcomp <- modelfit$aov.tab[-1,2]
        
        # Number of observations
        nvector <- nfun(modelfit)
        
        # Weights
        weights <- 1+(1/nvector)
        
        # weighted variance
        weighted_var_list <- Map("*", varcomp, weights)
        
        # Var(y) + Var(mean)
        pred_var <- sum(unlist(weighted_var_list))
        
        # Variances of variance components (without total)
        var_var <-modelfit$aov.tab$`Var(VC)`[-1]
        
        # Variance of variance components * weights^2 
        weighted_var_var <- Map("*", var_var, weights^2)
        
        # s^2(2), first formula in vdH pp.787
        pred_var_var <- sum(unlist(weighted_var_var))
        
        # dfs (according to vdH)
        dfs <-  max(1, min(nvector["error"]-1,
                           2*(pred_var^2/pred_var_var)))
        
        
        #-----------------------------------------------------------------------
        # Fixed effect (mu, intercept)
        
        muhat <- modelfit$FixedEffects@x
        
        #-----------------------------------------------------------------------
        # Calculation of the prediction bounds
        
        if(type=="upper"){
                
                upper <- muhat + qt(p=1-alpha, df=dfs)*sqrt(pred_var)   
                
                return(upper)
        }
        
        
        if(type=="lower"){
                
                lower <- muhat - qt(p=alpha, df=dfs)*sqrt(pred_var)   
                
                return(lower)
        }
        
        if(type=="both"){
                
                upper <- muhat + qt(p=1-alpha/2, df=dfs)*sqrt(pred_var)
                lower <- muhat - qt(p=1-alpha/2, df=dfs)*sqrt(pred_var)
                
                interval <- c(lower=lower, upper=upper)
                
                return(interval)
        }
        
}


#------------------------------------------------------------------------------

# Prediction interval using a model fit with remlMM from VCA
system.time(reml_pi_hb <- remlPI(modelfit=remlMM(log_nmru~1 + (run) + (mouse_id), dat_hb)))

# interval on the log scale
reml_pi_hb

# interval on response scale
exp(reml_pi_hb)







