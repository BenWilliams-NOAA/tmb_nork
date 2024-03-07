library(TMB)
library(tidyverse)
setwd(here::here('R'))
compile("nork.cpp")

# source data inputs from 2022 assessment
source(here::here('data', 'nork_data_2022.r'))
data = list(rec_age = rec_age,
            years = years,
            length_bins = length_bins,
            spawn_mo = spawn_mo,
            waa = waa,
            wt_mature = waa * maa,
            catch_obs = catch_obs,
            yr_catchwt = 1977,
            wt0 = 5,
            wt1 = 50,
            srv_ind = srv_ind,
            srv_obs = srv_obs,
            srv_sd = srv_sd,
            srv_wt = 0.25,
            fish_age_ind = fish_age_ind,
            fish_age_obs = fish_age_obs,
            fish_age_iss = fish_age_iss,
            fish_age_wt = 0.5,
            srv_age_ind = srv_age_ind,
            srv_age_obs = srv_age_obs,
            srv_age_iss = srv_age_iss,
            srv_age_wt = 0.5,
            fish_size_ind = fish_size_ind,
            fish_size_obs = fish_size_obs,
            fish_size_iss = fish_size_iss,
            fish_size_wt = 0.5,
            age_error = age_error,
            size_age = size_age,
            sigmaR = 1.5) # not currently used)
str(data) # check yo dims!!
pars = list(log_a50 = log_a50,
            delta = delta,
            log_mean_F = log_mean_F,
            log_Ft = log_Ft,
            log_M = log_M,
            log_mean_R = log_mean_R,
            log_Rt =  log_Rt,
            init_logRt = init_logRt,
            log_q = log_q)
 str(pars)           
dyn.load("nork")
obj = TMB::MakeADFun(data = data, 
                     parameters = pars,
                     DLL = "nork",
                     map = list(log_a50 = rep(factor(NA), 2),
                                delta = rep(factor(NA), 2),
                                # log_mean_F = factor(NA),
                                log_Ft = rep(factor(NA), T),
                                log_M = factor(NA),
                                log_mean_R = factor(NA),
                                log_Rt = rep(factor(NA), T),
                                init_logRt = rep(factor(NA), 48),
                                log_q = factor(NA)))
# Optimize the model
fit = nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr)
sd = sdreport(obj)
rep = obj$report()

rep$age_likelihood
rep$srv_age_likelihood

data.frame(years = srv_yrs,
obs = srv_obs,
pred = rep$srv_pred) %>% 
ggplot(aes(years, obs)) + 
geom_point() +
geom_line(aes(y=pred)) +
expand_limits(y = 0)



setwd(here::here())