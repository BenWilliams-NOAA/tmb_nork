library(TMB)
library(tidyverse)
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

try_improve = tryCatch(expr =
                         for(i in 1:2) {
                           g = as.numeric(obj$gr(fit$par))
                           h = optimHess(fit$par, fn = obj$fn, gr = obj$gr)
                           fit$par = fit$par - solve(h,g)
                           fit$objective = obj$fn(fit$par)
                         }
                       , error = function(e){e}, warning = function(w){w})

if(inherits(try_improve, "error") | inherits(try_improve, "warning")) {
  cat("didn't converge!!!\n")
  }

rep = obj$report(par = obj$env$last.par.best)

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


ages = rec_age:(nrow(fish_age_obs)+1)

data.frame(year = years,
          tot = rep$totbio,
          ssb = rep$ssb) %>% 
ggplot(aes(year, tot)) +
geom_point() + 
geom_line(aes(y=ssb)) +
expand_limits(y = 0) + 
geom_hline(yintercept = rep$Bzero * 0.4, lty = 3)


rep$fish_age_pred %>% 
        data.frame() %>% 
        mutate(age = ages) %>% 
        pivot_longer(-age) %>% 
        mutate(year = as.numeric(gsub("X", "", name)),
        type = 'pred')  %>% 
        bind_rows(fish_age_obs %>% 
        data.frame() %>% 
        mutate(age = ages) %>% 
        pivot_longer(-age) %>% 
        mutate(year = as.numeric(gsub("X", "", name)),
        type = 'obs')) %>% 
        ggplot(aes(x = age, y = value, color = type)) +
        geom_point() +
        facet_wrap(~year)

plot(years, rep$Nat[1,], type = 'bar')

 rep$srv_age_pred %>% 
        data.frame() %>% 
        mutate(age = ages) %>% 
        pivot_longer(-age) %>% 
        mutate(year = as.numeric(gsub("X", "", name)),
        type = 'pred')  %>% 
        bind_rows(srv_age_obs %>% 
        data.frame() %>% 
        mutate(age = ages) %>% 
        pivot_longer(-age) %>% 
        mutate(year = as.numeric(gsub("X", "", name)),
        type = 'obs')) %>% 
        ggplot(aes(age, value, color = type)) +
        geom_point() +
        facet_wrap(~year)   

rep$fish_size_pred %>% 
        data.frame() %>% 
        mutate(length = length_bins) %>% 
        pivot_longer(-length) %>% 
        mutate(year = as.numeric(gsub("X", "", name)),
        type = 'pred')  %>% 
        bind_rows(fish_size_obs %>% 
        data.frame() %>% 
        mutate(length = length_bins) %>% 
        pivot_longer(-length) %>% 
        mutate(year = as.numeric(gsub("X", "", name)),
        type = 'obs')) %>% 
        ggplot(aes(length, value, color = type)) +
        geom_point() +
        facet_wrap(~year, scales = 'free') 

