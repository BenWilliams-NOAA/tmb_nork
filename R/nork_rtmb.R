# attenpt at setting up northern rockfish model using RTMB
# load ----
# remotes::install_github("https://github.com/kaskr/RTMB", subdir="RTMB")
library(RTMB)
library(tidyverse)
library(Matrix)
library(tmbstan)
library(shinystan)
theme_set(afscassess::theme_report())

# data ----
# data and model
catch <- read.csv(here::here('data', 'catch.csv'))
srv <- read.csv(here::here('data', 'survey.csv'))
fac <- read.csv(here::here('data', 'fac.csv'))
sac <- read.csv(here::here('data', 'sac.csv'))
fsc <- read.csv(here::here('data', 'fsc.csv'))
slx <- read.csv(here::here('data', 'selex.csv'))
bio <- read.csv(here::here('data', 'bio_rec_f.csv'))
n_proj <- read.csv(here::here('data', 'n_proj.csv'))
b40 <- read.csv(here::here('data', 'b35_b40_yld.csv'))
source(here::here('R', 'model.r'))
source(here::here('data', 'nork_data_2022.r'))
source(here::here('R', 'utils.r'))
# inputs ------------------
data <- list(ages = ages,
             years = years,
             length_bins = 15:45,
             waa = waa,
             wt_mature = maa * waa / 2,
             spawn_mo = 5,
             catch_obs = catch_obs,
             catch_wt = c(rep(5, 17), rep(50, 45)),
             srv_obs = srv_obs,
             srv_ind = srv_ind,
             srv_sd = srv_sd,
             srv_wt = 0.25,
             fish_age_obs = fish_age_obs,
             fish_age_ind = fish_age_ind,
             fish_age_iss = fish_age_iss,
             fish_age_wt = 0.5,
             srv_age_obs = srv_age_obs,
             srv_age_ind = srv_age_ind,
             srv_age_iss = srv_age_iss,
             srv_age_wt = 0.5,
             fish_size_obs = fish_size_obs,
             fish_size_ind = fish_size_ind,
             fish_size_iss = fish_size_iss,
             fish_size_wt = 0.5,
             age_error = age_error,
             size_age = size_age,
             wt_fmort_reg = 0.1,
             wt_rec_var = 1,
             f_regularity = 0.1,
             mean_M = 0.06,
             cv_M = 0.05,
             mean_q = 1,
             cv_q = 0.45,
             mean_sigmaR = 1.5,
             cv_sigmaR = 0.01,
             yield_ratio = yield_ratio)
# parameters ----------------
pars <- list(log_M = log(0.05949983),
             log_a50C = log(8.23719578492),
             deltaC = 1.91866640868,
             log_a50S = log(9.09355629491),
             deltaS = 4.31923047888,
             log_q = log(0.8649374),
             log_mean_R = 3.50391,
             init_log_Rt = c(-1.31784331510, -1.41353004574, -1.45443795306, -1.44801600013, -1.44317223143,
                             -1.47227774217, -1.51491209802, -1.54201844034, -1.54900822105, -1.54357308125, -1.53271591530,
                             -1.51992650551, -1.50655074209, -1.49298490737, -1.47936419813, -1.46579792964, -1.45236085259,
                             -1.43911454732, -1.42609004023, -1.41334175243, -1.40089855851, -1.38878689192, -1.37702228101,
                             -1.36561634960, -1.35457991448, -1.34391631981, -1.33362760589, -1.32371409445, -1.31417494445,
                             -1.30500561365, -1.29620121156, -1.28775584994, -1.27966176612, -1.27190998138, -1.26449234040,
                             -1.25740132870, -1.25062646041, -1.24415861470, -1.23798650648, -1.23210097519, -1.22649216166,
                             -1.22115043981, -1.21606506031, -1.21122735382, -1.20662780636, -1.20225344500, -1.19809591673,
                             -1.19414750188),
             log_Rt = c(-1.18621056605,
                        -1.03273909462,-0.999797306673 ,-1.17237515921,-1.31481440189, -1.29947684085,
                        -1.19650952919, -1.07269785625, -0.858740244220 ,-0.648271337659, -0.687069961675,
                        0.445219166778, -0.468276695996, -0.563352814004, -0.453027729654, -0.584195538322,
                        -0.429755690388, 0.614844946205, 0.243664403047, -0.400119220417, -0.500707835112,
                        -0.335760494094, -0.00147160235074, 0.205169356072, -0.589731491602, 0.662087934070,
                        -0.0974724220631, -0.772529863247, -0.459399034102, -0.338932984917 ,-1.22764801352,
                        -0.470463642208, -0.821606342481, -0.909630584967, -1.23654078531, 0.436194651298,
                        0.0259279936320, -0.463748734096, -0.359703385166, 0.182589967526, -0.706622061830,
                        -1.01841380990, -0.793317388193, -1.45663999946, -2.24234182204, -2.25003316068,
                        -1.76689920898, -1.93202614806, -1.74725802806, -1.63462815593, -1.62908024719,
                        -1.84514613721, -2.16498687286, -1.76579171915, -1.89444424903, -1.60259969329,
                        -1.75681180416, -1.48983850684, -1.39173905744, -1.25892345656, -1.18461587351, -1.12500000434),
             log_mean_F = -3.58385236595,
             log_Ft =  c(-1.28306124683, 0.140199833981, 0.952641543234, 1.67793968633, 2.28974699704, 1.91823978977,
                         1.51557034959, 1.40878774713, 1.12294109172, 0.646529313912, 1.20801861960, 1.24452195275, 0.942561137836,
                         0.835801151139, 0.798114766724, 0.644693350645, -0.735539824209, -1.02062427087, -0.991743165230,
                         -0.920367489050, -0.403644808221, 0.516721388258, 0.389882559033, -0.981570055502, -2.79223891420,
                         -2.60309212227, -2.00728646156, -1.23024028545, -0.956209854061, -0.890429678442, 0.0349559024681,
                         0.546156698362, 0.0400591009128, 0.234941837688, 0.181355511206, -0.335989512042, -0.463551439451,
                         -0.416583500320, 0.167372898873, -0.293055078553, -0.341087697776, -0.279675928688, 0.155952440289,
                         0.0591776963077, -0.00563314517619, 0.0822175014391, -0.0885440790171, -0.113680406335, -0.120666350781,
                         -0.106971346707, -0.196222087322, 0.238469258157, 0.258217566983, 0.185811847296, 0.161828748432,
                         0.0782651773294, -0.500389551457, -0.207598679635, -0.00551620573419, -0.0956412655191, -0.049808327, -0.241028829),
             log_F35 = log(0.0736156723516),
             log_F40 = log(0.0612999005103),
             log_F50 = log(0.0428699169247),
             sigmaR = 1.5)

# map ----
map = list(log_M = factor(NA),
           log_a50C = factor(NA),
           deltaC = factor(NA),
           log_a50S = factor(NA),
           deltaS = factor(NA),
           log_q = factor(NA),
           log_mean_R = factor(NA),
           init_log_Rt = factor(rep(NA, 48)),
           log_Rt = factor(rep(NA, 62)),
           log_mean_F = factor(NA),
           log_Ft = factor(rep(NA, 62)),
           log_F35 = factor(NA),
           log_F40 = factor(NA),
           log_F50 = factor(NA),
           sigmaR = factor(NA))

# comparison run to admb model
obj_og <- RTMB::MakeADFun(f, 
                          pars, 
                          map = map)      
report_og <- obj_og$report(obj_og$env$last.par.best)
proj_og = proj_bio(report_og)

# compare to admb model
# look at the .rep file for more values for comparison
plot(years, bio$recruits)
lines(years, report_og$Nat[1,])

plot(years, bio$sp_biom)
lines(years, report_og$spawn_bio)

plot(years, bio$F)
lines(years, (report_og$Ft * max(report_og$slx[,1])))

report_og$B40
report_og$B35
b40

# run base model using output as start values, but with bounds
# sigmaR is not estimated in the GOA northern rockfish model
obj_1 <- RTMB::MakeADFun(f, 
                         pars, 
                         map = list(log_M = factor(NA),
                                    sigmaR = factor(NA)))  
length(obj_1$par)
names(obj_1$par)
lower = c(log(0.05), log(3), .5, log(3), 0.5, log(0.2), -10, rep(-30, length(pars$init_log_Rt)), 
          rep(-20, length(pars$log_Rt)), -20, rep(-30, length(years)), rep(-30,3))
upper = c(log(0.15), log(12), 5.5, log(12), 5.5, log(1.2), 10, rep(10, length(pars$init_log_Rt)), 
          rep(10, length(pars$log_Rt)), 10, rep(10, length(years)), rep(10,3))

fit_1 <- nlminb(obj_1$par,
                  obj_1$fn,
                  obj_1$gr,
                  control = list(iter.max=100000,
                                 eval.max=20000,
                                 abs.tol=0),
                  lower = lower,
                  upper = upper)

# has a hard time finding a global min
# $message
# [1] "relative convergence (4)"

sd_1 = RTMB::sdreport(obj_1)
# Hessian of fixed effects was not positive definite.
# some of the init_log_Rt don't have standard errors
report_1 <- obj_1$report(obj_1$env$last.par.best)
proj_bio(report_1)

# compare to admb model
plot(years, bio$recruits)
lines(years, report_1$Nat[1,])

plot(years, bio$sp_biom, ylim = c(10000, 85000))
lines(years, report_1$spawn_bio)

plot(years, bio$F, ylim = c(0, 0.4))
lines(years, (report_1$Ft * max(report_1$slx[,1])))

report_og$B40
report_1$B40

report_1$M
report_1$q



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
              sigmaR = 1.5)

obj_2 <- RTMB::MakeADFun(f, 
                         pars2, 
                         map = list(sigmaR = factor(NA)))  

fit_2 <- nlminb(obj_2$par,
                obj_2$fn,
                obj_2$gr,
                control = list(iter.max=100000,
                               eval.max=20000),
                lower = lower,
                upper = upper)
# has a hard time finding a global min
# Warning messages:
#   1: In nlminb(obj_2$par, obj_2$fn, obj_2$gr, control = list(iter.max = 1e+05,  :
#                                                                NA/NaN function evaluation
# gets to the same value as obj_1

sd_2 = RTMB::sdreport(obj_2)
# similar results as previous run
report_2 <- obj_2$report(obj_2$env$last.par.best)
proj_bio(report_2)

# compare to admb model
plot(years, bio$recruits)
lines(years, report_1$Nat[1,])
lines(years, report_2$Nat[1,], col=3)

plot(years, bio$sp_biom, ylim = c(10000, 85000))
lines(years, report_1$spawn_bio)
lines(years, report_2$spawn_bio, col=3)

plot(years, bio$tot_biom, ylim = c(30000, 200000))
lines(years, report_1$tot_bio)
lines(years, report_2$tot_bio, col=3)

plot(years, bio$F, ylim = c(0, 0.4))
lines(years, (report_1$Ft * max(report_1$slx[,1])))
lines(years, (report_2$Ft * max(report_2$slx[,1])), col = 3)

report_og$B40
report_1$B40
report_2$B40

report_1$Nat - report_2$Nat
report_1$spawn_bio - report_2$spawn_bio
report_1$tot_bio - report_2$tot_bio
report_1$Ft - report_2$Ft


# don't really like older catch values as they are just expansions from pop catches
# remove anything before 1979 - to match what we do with duskies

years = 1979:2022
catch_obs = c(670, 810, 1477,
              3920, 3618, 1002, 185, 248, 483, 1107, 1527, 1716, 4528, 7770 ,4820.3,
              5966.19, 5635.4, 3340.31, 2935.01, 3055.45, 5408.89, 3333.47, 3132.77, 3339.25, 5256.072,
              4810.786, 4521.671, 4957.706, 4186.649, 4052.09, 3951.653, 3902.3631, 3443.9287, 5077.0064, 4879.2848, 4278.2316, 3944.5536,
              3433.9938, 1835.0952, 2358.7904, 2748.1593, 2384.9172, 2376.4923, 1875.7379)

srv_ind = ifelse(years %in% srv_yrs, 1, 0)
fish_age_ind = ifelse(years %in% fish_age_yrs, 1, 0)
srv_age_ind = ifelse(years %in% srv_age_yrs, 1, 0)
fish_size_ind = ifelse(years %in% fish_size_yrs, 1, 0)

data <- list(ages = ages,
             years = years,
             length_bins = 15:45,
             waa = waa,
             wt_mature = maa * waa / 2,
             spawn_mo = 5,
             catch_obs = catch_obs,
             catch_wt = rep(50, length(years)),
             srv_obs = srv_obs,
             srv_ind = srv_ind,
             srv_sd = srv_sd,
             srv_wt = 0.25,
             fish_age_obs = fish_age_obs,
             fish_age_ind = fish_age_ind,
             fish_age_iss = fish_age_iss,
             fish_age_wt = 0.5,
             srv_age_obs = srv_age_obs,
             srv_age_ind = srv_age_ind,
             srv_age_iss = srv_age_iss,
             srv_age_wt = 0.5,
             fish_size_obs = fish_size_obs,
             fish_size_ind = fish_size_ind,
             fish_size_iss = fish_size_iss,
             fish_size_wt = 0.5,
             age_error = age_error,
             size_age = size_age,
             wt_fmort_reg = 0.1,
             mean_M = 0.06,
             sd_M = 0.05,
             mean_q = 1,
             sd_q = 0.45,
             mean_a50C = log(8),
             sd_a50C = 0.25,
             mean_deltaC = 3,
             sd_deltaC = 1,
             mean_a50S = log(8),
             sd_a50S = 0.25,
             mean_deltaS = 3,
             sd_deltaS = 1,
             mean_sigmaR = 1.5,
             sd_sigmaR = 0.01,
             yield_ratio = yield_ratio)

# run base using orig input values 
pars3 <- list(log_M = log(0.0595),
              log_a50C = log(7.5),
              deltaC = 3,
              log_a50S = log(7.3),
              deltaS = 3.8,
              log_q = 0,
              log_mean_R = 4.3,
              init_log_Rt = rep(0, 48),
              log_Rt = rep(0, length(years)),
              log_mean_F = 0,
              log_Ft =  rep(0, length(years)),
              log_F35 = 0,
              log_F40 = 0,
              log_F50 = 0,
              sigmaR = 1.5)

obj_3 <- RTMB::MakeADFun(f, 
                         pars3, 
                         map = list(sigmaR = factor(NA)))  

lower = c(log(0.05), log(3), .5, log(3), 0.5, log(0.2), -15, rep(-15, length(pars3$init_log_Rt)), 
          rep(-15, length(pars3$log_Rt)), -15, rep(-15, length(years)), rep(-15,3))
upper = c(log(0.15), log(12), 5.5, log(12), 5.5, log(1.2), 10, rep(15, length(pars3$init_log_Rt)), 
          rep(15, length(pars3$log_Rt)), 10, rep(10, length(years)), rep(10,3))

fit_3 <- nlminb(obj_3$par,
                obj_3$fn,
                obj_3$gr,
                control = list(iter.max=100000,
                               eval.max=20000),
                lower = lower,
                upper = upper)

sd_3 = RTMB::sdreport(obj_3)
# similar results as previous run
report_3 <- obj_3$report(obj_3$env$last.par.best)
proj_bio(report_2)
proj_bio(report_3)

# compare to admb model
plot(years, report_3$Nat[1,], type ='l')
plot(years, report_3$spawn_bio, type ='l')
plot(years, report_3$tot_bio, type ='l')
lines(1961:2022, bio$tot_biom, ylim = c(30000, 200000))
lines(1961:2022, report_1$tot_bio)

plot(1961:2022, bio$F, ylim = c(0, 0.4))
lines(1961:2022, (report_1$Ft * max(report_1$slx[,1])))
lines(years, (report_3$Ft * max(report_3$slx[,1])), col = 3)

report_og$M
report_1$M
report_2$M
report_3$M

report_og$q
report_1$q
report_2$q
report_3$q
# q is quite low, tighten up the prior

data <- list(ages = ages,
             years = years,
             length_bins = 15:45,
             waa = waa,
             wt_mature = maa * waa / 2,
             spawn_mo = 5,
             catch_obs = catch_obs,
             catch_wt = rep(50, length(years)),
             srv_obs = srv_obs,
             srv_ind = srv_ind,
             srv_sd = srv_sd,
             srv_wt = 0.25,
             fish_age_obs = fish_age_obs,
             fish_age_ind = fish_age_ind,
             fish_age_iss = fish_age_iss,
             fish_age_wt = 0.5,
             srv_age_obs = srv_age_obs,
             srv_age_ind = srv_age_ind,
             srv_age_iss = srv_age_iss,
             srv_age_wt = 0.5,
             fish_size_obs = fish_size_obs,
             fish_size_ind = fish_size_ind,
             fish_size_iss = fish_size_iss,
             fish_size_wt = 0.5,
             age_error = age_error,
             size_age = size_age,
             wt_fmort_reg = 0.1,
             mean_M = 0.06,
             sd_M = 0.05,
             mean_q = 1,
             sd_q = 0.15,
             mean_sigmaR = 1.5,
             sd_sigmaR = 0.01,
             yield_ratio = yield_ratio)


obj_3 <- RTMB::MakeADFun(f, 
                         pars3, 
                         map = list(sigmaR = factor(NA)))  

fit_3 <- nlminb(obj_3$par,
                obj_3$fn,
                obj_3$gr,
                control = list(iter.max=100000,
                               eval.max=20000),
                lower = lower,
                upper = upper)

sd_3 = RTMB::sdreport(obj_3)
# similar results as previous run
report_3 <- obj_3$report(obj_3$env$last.par.best)
proj_bio(report_og)
proj_bio(report_2)
proj_bio(report_3)

report_2$q
report_3$q

report_2$M
report_3$M

report_og$B40
report_2$B40
report_3$B40

report_og$a50C
report_2$a50C
report_3$a50C

report_og$a50S
report_2$a50S
report_3$a50S

plot(1961:2022, report_og$Nat[1,], ylim = c(0,67), type = 'l', col='blue')
lines(1961:2022, report_2$Nat[1,], ylim = c(0,67), type = 'l', col='darkgray')
lines(years, report_3$Nat[1,])

plot(1961:2022, report_2$spawn_bio, ylim = c(10000, 85000), type = 'l', col='darkgray')
lines(years, report_3$spawn_bio)

plot(1961:2022, report_2$tot_bio, ylim = c(10000, 210000), type = 'l', col='darkgray')
lines(years, report_3$tot_bio)

plot(srv_yrs, report_2$srv_pred, type = 'l', ylim = c(60000, 150000), col='darkgray')
lines(srv_yrs, report_3$srv_pred)

plot(1961:2022, (report_2$Ft * max(report_1$slx[,1])), type = 'l', col='darkgray')
lines(years, (report_3$Ft * max(report_1$slx[,1])))

# mcmc ----
# not gonna work out with this limited number of chaines
mcmc <- tmbstan(obj_3, chains=1, iter=4000)
post <- as.matrix(mcmc)
# don't run unless you have some time on your hands...
sims <- proj_bio(report_3, obj_3, post, reps=nrow(post))

ggplot(sims, aes(year, value, color = id, group = id)) +
  stat_summary(fun.y=mean, geom='line') +
  stat_summary(fun.data = mean_cl_boot, geom = "ribbon", alpha = 0.2, color = NA) +
  expand_limits(y = 0) +
  geom_hline(yintercept = c(report_3$B35, report_3$B40))

sb = get_bio(post, obj_3, reps = nrow(post))
post$spawn_bio
obj_3$report(post[2000,-ncol(post)])$spawn_bio
sb %>% 
  as.data.frame() %>% 
  mutate(year = report_3$years) %>% 
  tidyr::pivot_longer(-year) %>% 
  mutate(sim = as.numeric(gsub("V", "", name)),
         id = 'spawn_bio') %>% 
  dplyr::select(-name) %>% 
  bind_rows(sims) %>%
  dplyr::filter(id=='spawn_bio') %>% 
  ggplot(aes(year, value)) +
  stat_summary(fun.y=mean, geom='line') +
  stat_summary(fun.data = mean_cl_normal, geom = "ribbon", alpha = 0.2, color = NA)+
  expand_limits(y = 0) +
  geom_hline(yintercept = c(report_3$B35, report_3$B40))

# Choose parameter values to profile
par_names <- c("log_a50C", "log_a50S")  # Replace with your parameter names
par_ranges <- list(
  log_a50C = seq(-3, 3, length.out=50),  # Replace with appropriate range
  log_a50S = seq(-3, 3, length.out=50)   # Replace with appropriate range
)

# Profile the likelihood
profile_results <- multi_likelihood(par_names, par_ranges, obj_3, fit_3)
ggplot(profile_results, aes(log_a50C, log_a50S)) +
  geom_tile(aes(fill = log_like)) +
  stat_contour(aes(z = log_like)) +
  scico::scale_fill_scico(direction = -1, palette = 'oslo')
scale_fill_brewer(
  type = "seq",
  palette = "blues",
  direction = -1,
)

par_name = 'log_a50C'
par_values = seq(-3, 3, length.out = 50)
par_like = single_likelihood(par_name, par_values, obj_3, fit_3)
# Plot the profile
plot(par_values, par_like, type="l", xlab=par_name, ylab="Log-Likelihood")

par_name = 'log_a50S'
par_values = seq(0, 3, length.out = 50)
par_like = single_likelihood(par_name, par_values, obj_3, fit_3)
# Plot the profile
plot(par_values, par_like, type="l", xlab=par_name, ylab="Log-Likelihood")

par_name = 'log_mean_R'
par_values = seq(2.2, 2.6, length.out = 50)
par_like = single_likelihood(par_name, par_values, obj_3, fit_3)
# Plot the profile
plot(par_values, par_like, type="l", xlab=par_name, ylab="Log-Likelihood")
