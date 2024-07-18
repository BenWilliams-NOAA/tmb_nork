# attenpt at setting up northern rockfish model using RTMB
# remotes::install_github("https://github.com/kaskr/RTMB", subdir="RTMB")
library(RTMB)
library(tidyverse)
library(Matrix)
library(tmbstan)
library(shinystan)
theme_set(afscassess::theme_report())
# data and model
source('orig_pars.r')
catch <- read.csv(here::here('data', 'catch.csv'))
srv <- read.csv(here::here('data', 'survey.csv'))
fac <- read.csv(here::here('data', 'fac.csv'))
sac <- read.csv(here::here('data', 'sac.csv'))
fsc <- read.csv(here::here('data', 'fsc.csv'))
slx <- read.csv(here::here('data', 'selex.csv'))
bio <- read.csv(here::here('data', 'bio_rec_f.csv'))
n_proj <- read.csv(here::here('data', 'n_proj.csv'))
b40 <- read.csv(here::here('data', 'b35_b40_yld.csv'))
source('model.r')

# run base using orig pars as input values 
obj5 <- RTMB::MakeADFun(f, 
                        pars, 
                        map = list(sigmaR = factor(NA)))   
fit5 <- nlminb(obj5$par,
               obj5$fn,
               obj5$gr)
fit5 = TMBhelper::fit_tmb(obj5,
                          newtonsteps = 1,
                          get_sd = FALSE)
sd5 = RTMB::sdreport(obj5)
report5 <- obj5$report(obj5$env$last.par.best)

proj_bio(report5)
report5$M
report5$q
report5$B40

par_name = 'log_mean_R'
par_values = seq(1.6, 3, length.out = 50)
par_like = single_likelihood(par_name, par_values, obj5, fit5)
plot(par_values, par_like, type="l", xlab=par_name, ylab="Log-Likelihood")

par_name = 'log_M'
par_values = seq(-10, 0, length.out = 50)
par_like = single_likelihood(par_name, par_values, obj5, fit5)
plot(par_values, par_like, type="l", xlab=par_name, ylab="Log-Likelihood")

par_name = 'log_q'
par_values = seq(-4, 4, length.out = 50)
par_like = single_likelihood(par_name, par_values, obj5, fit5)
plot(par_values, par_like, type="l", xlab=par_name, ylab="Log-Likelihood")

# run base using orig pars as input values  - but with updated prior on M
data$mean_M <- 0.07
data$sd_M <- 0.000873

data$sd_q <- 0.45

obj6 <- RTMB::MakeADFun(f, 
                        pars, 
                        map = list(sigmaR = factor(NA),
                                   log_M = factor(NA)))   
fit6 <- nlminb(obj6$par,
               obj6$fn,
               obj6$gr)
# sd6 = RTMB::sdreport(obj6)
report6 <- obj6$report(obj6$env$last.par.best)
proj_bio(report6)

report6$M
report6$q
report6$B40

par_name = 'log_mean_R'
par_values = seq(1.6, 3, length.out = 50)
par_like = single_likelihood(par_name, par_values, obj5, fit5)
plot(par_values, par_like, type="l", xlab=par_name, ylab="Log-Likelihood")

par_name = 'log_q'
par_values = seq(-4, 4, length.out = 50)
par_like = single_likelihood(par_name, par_values, obj6, fit6)
plot(par_values, par_like, type="l", xlab=par_name, ylab="Log-Likelihood")


###################
data$mean_M = 0.06
data$sd_M = 0.05

# run base using orig input values 
pars2 <- list(log_M = log(0.0595),
              log_a50C = log(7.5),
              deltaC = 3,
              log_a50S = log(7.3),
              deltaS = 3.8,
              log_q = 0,
              log_mean_R = 4.3,
              init_log_Rt = rep(0, 48),
              log_Rt = rep(0, 62),
              log_mean_F = 0,
              log_Ft =  rep(0, 62),
              log_F35 = 0,
              log_F40 = 0,
              log_F50 = 0,
              sigmaR = 1.5,
              dum = 0)
obj7 <- RTMB::MakeADFun(f, 
                        pars2, 
                        map = list(log_M = factor(NA),
                                   sigmaR = factor(NA)))   
fit7 <- nlminb(obj7$par,
               obj7$fn,
               obj7$gr)
# fit7 = TMBhelper::fit_tmb(obj7,
#                           newtonsteps = 1,
#                           get_sd = FALSE)
sd7 = RTMB::sdreport(obj7)
report7 <- obj7$report(obj7$env$last.par.best)
proj_bio(report7)
report7$M
report7$q




# selectivity block ----
pars3 <- list(log_M = log(0.0595),
              log_a50C = log(7.5),
              deltaC = 3,
              log_a50S = log(7.3),
              deltaS = 3.8,
              log_q = 0,
              log_mean_R = 4.3,
              init_log_Rt = rep(0, 48),
              log_Rt = rep(0, 44),
              log_mean_F = 0,
              log_Ft =  rep(0, 44),
              log_F35 = 0,
              log_F40 = 0,
              log_F50 = 0,
              sigmaR = 1.5,
              dum = 0)
data$slx_blk = ifelse(years<=1997, 1, 2)
pars3$log_a50C2 = 3
pars3$deltaC2 = 3
data$mean_M
data$sd_M = 0.005
data$sd_q = 0.15
obj_blk <- RTMB::MakeADFun(f1, 
                        pars3, 
                        map = list(sigmaR = factor(NA)))   
fit_blk <- nlminb(obj_blk$par,
                  obj_blk$fn,
                  obj_blk$gr)
fit_blk = TMBhelper::fit_tmb(obj_blk,
                          newtonsteps = 1,
                          get_sd = FALSE)
report_blk <- obj_blk$report(obj_blk$env$last.par.best)
proj_bio(report_blk)
report_blk$M
report_blk$q
report_blk$slx
report_blk$catch_pred


### Fishery size compositions

data.frame(length = length_bins, 
           report_blk$fish_size_pred) %>% 
  pivot_longer(-length) %>% 
  mutate(year = rep(fish_size_yrs, length(length_bins)),
         groups = 'pred2') %>% 
  bind_rows(fsc) -> df4

df4 %>% 
  ggplot(aes(length, value, color = groups)) + 
  geom_line() +
  facet_wrap(~year)


df4 %>% 
  select(length, value, year, groups) %>% 
  filter(groups!='obs') %>% 
  pivot_wider(names_from=groups, values_from = value) %>% 
  mutate(diff = (pred2 - pred) / pred,
         Length = factor(length)) %>% 
  pivot_longer(diff) %>% 
  ggplot(aes(year, value, color = Length)) + 
  geom_line() +
  # facet_wrap(~year) +
  scale_y_continuous(labels = scales::percent) +
  ylab('percent difference') +
  geom_hline(yintercept=0, lty=3)