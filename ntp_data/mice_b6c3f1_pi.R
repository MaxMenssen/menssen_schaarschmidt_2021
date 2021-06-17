#------------------------------------------------------------------------------
#---------------- PI for m=2 based on HCD -------------------------------------
#------------------------------------------------------------------------------

# Data import
dat_b <- read.csv("ntp_data/mice_b6c3f1_female_battelle.csv")
dat_is <- read.csv("ntp_data/mice_b6c3f1_female_iit_southern.csv")

#------------------------------------------------------------------------------

# h2 model with lme4

library(lme4)

fit_b <- lmer(mmwbw~1 + (1|lab_name) + (1|lab_name:pathway), data=dat_b)
summary(fit_b)



#------------------------------------------------------------------------------

# prediction interval m=2

# Please note, that the columns of newdat have to be exactly the same as 
fit_b@frame

# Hence dat_is is reordered with
library(tidyverse)
select(dat_is,
       mmwbw,
       lab_name,
       pathway)


library(predint)
pi_m2 <- lmer_pi(model=fit_b, 
                 newdat=select(dat_is,
                               mmwbw,
                               lab_name,
                               pathway))

pi_m2

#   mmwbw                    lab_name         pathway hist_mean quant_calib  pred_se    lower    upper cover
# 1  62.6      IIT Research Institute         wbe_air  59.68749    2.663594 5.879907 44.02581 75.34918  TRUE
# 2  57.7 Southern Research Institute gavage_corn oil  59.68749    2.663594 5.879907 44.02581 75.34918  TRUE

# Since the PI is based on bootstrap-calibration, the output may change from
# run to run.


