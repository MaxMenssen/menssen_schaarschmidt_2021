#------------------------------------------------------------------------------
#------------------- Bootstrap-calibrated PI ----------------------------------
#------------------------------------------------------------------------------

# Data of Table 4 from Hoffmann & Berger 2011 (Normalized mean RU)

run_1 <- c(1.11, 1.00, 1.50, 1.04, 0.929, 1.02, 0.965, 0.973, 1.01, 1.02,
           1.04, 1.13, 1.03, 1.04, 1.02, 1.04, 0.929, 1.04, 0.947, 0.947)

run_2 <- c(1.25, 1.18, 1.52, 1.15, 1.13, 1.03, 1.18, 1.18, 1.21, 1.29,
           1.16, 1.22, 1.35, 1.28, 1.29, 1.38, 1.20, 1.16, 1.20, 1.10)

run_3 <- c(1.03, 0.961, 1.30, 1.12, 1.06, 1.10, 0.804, 0.850, 1.02, 1.07, 
           0.933, 0.979, 1.11, 1.12, 1.02, 1.03, 0.915, 1.02, 0.915, 0.961)

run <- factor(c(rep(1, 20), rep(2, 20), rep(3, 20)))

mouse_id <- factor(c(1:20, 1:20, 1:20))

# Data is given on log-scale
dat_hb <- data.frame(run=run,
                     mouse_id=mouse_id,
                     log_nmru=log(c(run_1, run_2, run_3)))

#------------------------------------------------------------------------------

library(lme4)
library(predint)

# random effects model with lmer
fit_hb_lmer <- lmer(log_nmru~1 + (1|run) + (1|mouse_id), data=dat_hb)

# prediction interval for one future observations
system.time(calib_pi <- lmer_pi(model=fit_hb_lmer,
                                m=1))

calib_pi



