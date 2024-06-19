# attenpt at setting up northern rockfish model using RTMB
install.packages('RTMB', repos = c('https://kaskr.r-universe.dev', 'https://cloud.r-project.org'))
library(Matrix)
library(tmbstan)
library(shinystan)

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
             mean_M = 0.07,
             sd_M = 0.015,
             mean_q = 1,
             sd_q = 0.45,
             mean_sigmaR = 1.5,
             sd_sigmaR = 0.1,
             proj_rec_yrs = proj_rec_yrs,
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
             sigmaR = 1.5,
             dum = 0)

# model ----
f <- function(pars) {
  RTMB::getAll(pars, data)
  # setup -------------
  M = exp(log_M)
  a50C = exp(log_a50C)
  a50S = exp(log_a50S)
  q = exp(log_q)
  spawn_fract = (spawn_mo - 1) / 12
  spawn_adj = exp(-M)^(spawn_fract)
  A = nrow(age_error)
  A1 = length(ages)
  T = length(catch_obs)
  Ts = sum(srv_ind)
  Tfa = sum(fish_age_ind)
  Tsa = sum(srv_age_ind)
  Tfs = sum(fish_size_ind)
  Tproj = 15
  L = length(length_bins)
  cwt = unique(catch_wt) # time blocks for catch weighting
  g = 0.00001 # small number
  F50 = exp(log_F50)
  F40 = exp(log_F40)
  F35 = exp(log_F35)
  f_regularity = 0.1
  
  # containers ---------
  slx = matrix(NA, A, 2)
  Bat = Cat = Nat = Fat = Zat = Sat = matrix(0, A, T)
  initNat = rep(0, A)
  catch_pred = rep(0, T)
  srv_pred = rep(0, Ts)
  fish_age_pred = matrix(0, A1, Tfa)
  srv_age_pred = matrix(0, A1, Tsa)
  fish_size_pred = matrix(0, L, Tfs)
  spawn_bio = tot_bio = rep(0, T)
  N_spr = sb_spr = matrix(1, A, 4)
  # N_proj = Cat_proj = matrix(0, A, Tproj)
  # spawn_bio_proj = tot_bio_proj = rep(0, Tproj)
  
  # priors -----------------
  nll_M = nll_q = 0.0
  nll_M = dnorm(M, mean_M, sd_M, TRUE)
  nll_q = dnorm(q, mean_q, sd_q, TRUE)
  nll_sigmaR = dnorm(sigmaR, mean_sigmaR, sd_sigmaR, TRUE)
  
  # analysis ----------
  # selectivity 
  sel <- function(a, a50, delta) {
    1. / (1. + exp(-2.944438979 * ((a+1) - a50) / delta))
  }
  
  for(a in 1:A) {
    slx[a,1] = sel(a, a50C, deltaC)
    slx[a,2] = sel(a, a50S, deltaS)
  }
  
  # mortality 
  Ft = exp(log_mean_F + log_Ft)
  for(t in 1:T){
    for(a in 1:A) {
      Fat[a,t] = Ft[t] * slx[a,1]
      Zat[a,t] = Fat[a,t] + M
      Sat[a,t] = exp(-Zat[a,t])
    }
  }
  
  # population ----
  ## Bzero
  # initNat[1] = exp(log_mean_R)
  # for (a in 2:A) {
  #   initNat[a] <- initNat[a-1] * exp(-M)
  # }   
  # initNat[A] <- initNat[A] / (1 - exp(-M))
  # Bzero = sum(initNat * wt_mature * spawn_adj)
  
  ## numbers at age
  # populate first row
  # need to use correct log_Rt to match ADMB model
  for(t in 1:T) {
    Nat[1,t] = exp(log_mean_R + log_Rt[t])
  }
  # populate first column
  for(a in 2:(A-1)) {
    Nat[a,1] = exp(log_mean_R - (a-1) * M + init_log_Rt[a-1])
  }
  Nat[A,1] = exp(log_mean_R - (A-1) * M) / (1 - exp(-M))
  
  for(t in 2:T) {
    for(a in 2:A) {
      Nat[a,t] = Nat[a-1,t-1] * Sat[a-1,t-1]
    }
    Nat[A,t] = Nat[A,t] + Nat[A,t-1] * Sat[A,t-1]  
  }
  
  #     # spawn_bio
  for(t in 1:T) {
    spawn_bio[t] = sum(Nat[,t] * wt_mature)
    tot_bio[t] = sum(Nat[,t] * waa)
  }
  # flag - should be:
  # spawn_bio[T] = sum(Nat[,T] * spawn_adj * wt_mature)
  
  
  # catch
  for(t in 1:T){
    for(a in 1:A){
      Cat[a,t] = Fat[a,t] / Zat[a,t] * Nat[a,t] * (1.0 - Sat[a,t])
    }
    catch_pred[t] = sum(Cat[,t] * waa)
  }
  
  # survey biomass
  isrv = 1
  srv_like = 0.0
  # survey biomass & likelihood
  for(t in 1:T) {
    if(srv_ind[t]==1) {
      # for(a in 1:A) {
      srv_pred[isrv] = sum(Nat[,t] * slx[,2] * waa)
      # }
      srv_pred[isrv] = srv_pred[isrv] * q
      srv_like = srv_like + sum((log(srv_obs[isrv]) - log(srv_pred[isrv]))^2 /
                                  (2 * (srv_sd[isrv] / srv_obs[isrv])^2))
      isrv = isrv + 1
    }
  }
  
  # fishery age comp
  icomp = 1
  for(t in 1:T) {
    if(fish_age_ind[t] == 1) {
      fish_age_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * age_error))
      # fish_age_lk[icomp] = sum(fish_age_iss[icomp] * ((fish_age_obs[,icomp] + g) * log((fish_age_pred[,icomp] + g / (fish_age_obs[,icomp] + g)))))
      icomp = icomp + 1
    }
  }
  
  # survey age comp
  icomp = 1
  for(t in 1:T) {
    if(srv_age_ind[t] == 1) {
      srv_age_pred[,icomp] = as.numeric(colSums((Nat[,t] * slx[,2]) / sum(Nat[,t] * slx[,2]) * age_error))
      icomp = icomp + 1
    }
  }
  
  # fishery size comp
  icomp = 1
  for(t in 1:T) {
    if(fish_size_ind[t] == 1) {
      fish_size_pred[,icomp] = as.numeric(colSums((Cat[,t] / sum(Cat[,t])) * size_age))
      # fish_age_lk[icomp] = sum(fish_age_iss[icomp] * ((fish_age_obs[,icomp] + g) * log((fish_age_pred[,icomp] + g / (fish_age_obs[,icomp] + g)))))
      icomp = icomp + 1
    }
  }
  
  # SPR ------------------------
  data.frame(ind = proj_rec_yrs,
             log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # filter years 1979:
  df = df[1:(T-ages[1]),]
  df = df[df$ind>=1,]
  pred_rec = mean(df$pred_rec)
  
  for(a in 2:A) {
    N_spr[a,1] = N_spr[a-1,1] * exp(-M)
    N_spr[a,2] = N_spr[a-1,2] * exp(-(M + F50 * slx[a-1,1]))
    N_spr[a,3] = N_spr[a-1,3] * exp(-(M + F40 * slx[a-1,1]))
    N_spr[a,4] = N_spr[a-1,4] * exp(-(M + F35 * slx[a-1,1]))
  }
  
  N_spr[A,1] = N_spr[A-1,1] * exp(-M) / (1 - exp(-M))
  N_spr[A,2] = N_spr[A-1,2] * exp(-(M + F50 * slx[A-1,1])) / (1 - exp(-(M + F50 * slx[A,1])))
  N_spr[A,3] = N_spr[A-1,3] * exp(-(M + F40 * slx[A-1,1])) / (1 - exp(-(M + F40 * slx[A,1])))
  N_spr[A,4] = N_spr[A-1,4] * exp(-(M + F35 * slx[A-1,1])) / (1 - exp(-(M + F35 * slx[A,1])))
  
  for(a in 1:A) {
    sb_spr[a,1] = N_spr[a,1] * wt_mature[a] * exp(-spawn_fract * M)
    sb_spr[a,2] = N_spr[a,2] * wt_mature[a] * exp(-spawn_fract * (M + F50 * slx[a,1]))
    sb_spr[a,3] = N_spr[a,3] * wt_mature[a] * exp(-spawn_fract * (M + F40 * slx[a,1]))
    sb_spr[a,4] = N_spr[a,4] * wt_mature[a] * exp(-spawn_fract * (M + F35 * slx[a,1]))
  }
  
  SB0 = sum(sb_spr[,1])
  SBF50 = sum(sb_spr[,2])
  SBF40 = sum(sb_spr[,3])
  SBF35 = sum(sb_spr[,4])
  
  sprpen = 100. * (SBF50 / SB0 - 0.5)^2
  sprpen = sprpen + 100. * (SBF40 / SB0 - 0.4)^2
  sprpen = sprpen + 100. * (SBF35 / SB0 - 0.35)^2
  
  B0 = SB0 * pred_rec
  B40 = SBF40 * pred_rec
  B35 = SBF35 * pred_rec
  
  # projection ----------------------
  # storage
  N_proj = Cat_proj = Cat_ofl_proj = Zabc_proj = Zofl_proj = S_proj = matrix(0, A, Tproj)
  catch_proj = catch_ofl_proj = spawn_bio_proj = tot_bio_proj = catch_proj = rep(0, Tproj)
  
  F40_proj = F40
  F35_proj = F35
  
  # total F
  Fabc_tot_proj = slx[,1] * F40_proj
  Fofl_tot_proj = slx[,1] * F35_proj
  
  # every years starts with mean recruits from 1979:(T-ages[1])
  N_proj[1,] = pred_rec
  
  # first year of proj
  for(a in 2:A) {
    N_proj[a,1] = Nat[a-1,T] * Sat[a-1,T]
  }
  N_proj[A,1] = Nat[A-1,T] * Sat[A-1,T] + Nat[A,T] * Sat[A,T]
  spawn_bio_proj[1] = sum(N_proj[,1] * exp(-yield_ratio * Fabc_tot_proj - M)^spawn_fract * wt_mature)
  
  # tier check
  # if((spawn_bio_proj[1] / B40) > 1) {
  #   F40_proj = F40
  #   F35_proj = F35
  #  } else {
  #   F40_proj[,1] = F40_proj * (spawn_bio_proj[1] / B402 - 0.05) / 0.95
  #   F35_proj[,1] = F35_proj * (spawn_bio_proj[1] / B402 - 0.05) / 0.95
  #  }
  Fabc_tot_proj = F40_proj * slx[,1]
  Fofl_tot_proj = F35_proj * slx[,1]
  Zabc_proj = Fabc_tot_proj + M
  Zofl_proj = Fofl_tot_proj + M
  S_proj = exp(-Zabc_proj)
  
  for(t in 2:Tproj) {
    for(a in 2:A){
      N_proj[a,t] = N_proj[a-1,t-1] * exp(-yield_ratio * Fabc_tot_proj[a] - M)
    }
    N_proj[A,t] = N_proj[A-1,t-1] * exp(-yield_ratio * Fabc_tot_proj[a] - M) +
      N_proj[A,t-1] * exp(-yield_ratio * Fabc_tot_proj[a] - M)
    
    spawn_bio_proj[t] = sum(N_proj[,t] * exp(-yield_ratio * Fabc_tot_proj[a] - M)^spawn_fract * wt_mature)
  }
  
  #   # tier check
  #   if(spawn_bio_proj[t] / B40 > 1) {
  #     F40_proj = F40
  #     F35_proj = F35
  #   } else {
  #     Fabc_proj[,t] = F40_proj * (spawn_bio_proj[t] / B40 - 0.05) / 0.95 
  #     Fofl_proj[,t] = F35_proj * (spawn_bio_proj[t] / B40 - 0.05) / 0.95 
  #   }
  #   Fabc_tot_proj[,t] = F40_proj * slx[,1]
  #   Fofl_tot_proj[,t] = F35_proj * slx[,1]
  #   Zabc_proj[,t] = Fabc_tot_proj[,t] + M
  #   Zofl_proj[,t] = Fofl_tot_proj[,t] + M
  #   S_proj[,t] = exp(-Zabc_proj[,t])
  #   
  # }
  # 
  for(t in 1:Tproj) {
    Cat_proj[,t] = yield_ratio * N_proj[,t] * Fabc_tot_proj / Zabc_proj * (1 - S_proj)
    Cat_ofl_proj[,t] = yield_ratio * N_proj[,t] * Fofl_tot_proj / Zofl_proj * (1 - exp(-Zofl_proj))
  }
  
  catch_proj = colSums(Cat_proj * waa / yield_ratio)
  catch_ofl_proj = colSums(Cat_ofl_proj * waa / yield_ratio)
  tot_bio_proj = colSums(N_proj * waa)
  
  # likelihoods --------------------
  # catch
  ssqcatch = sum(catch_wt * (log(catch_obs + g) - log(catch_pred + g))^2)
  
  # fishery age comp
  fish_age_lk = 0.0
  offset = 0.0
  for(t in 1:Tfa) {
    offset = offset - fish_age_iss[t] * sum((fish_age_obs[,t] + g) * log(fish_age_obs[,t] + g))
    fish_age_lk = fish_age_lk - sum(fish_age_iss[t] * (fish_age_obs[,t] + g) * log(fish_age_pred[,t] + g))
  }
  fish_age_lk = fish_age_lk - offset
  
  # survey age comp
  srv_age_lk = 0.0
  offset_sa = 0.0
  for(t in 1:Tsa) {
    # for(a in 1:A1) {
    offset_sa = offset_sa - srv_age_iss[t] * sum((srv_age_obs[,t] + g) * log(srv_age_obs[,t] + g))
    srv_age_lk = srv_age_lk - srv_age_iss[t] * sum((srv_age_obs[,t] + g) * log(srv_age_pred[,t] + g))
    # }
  }
  srv_age_lk = srv_age_lk - offset_sa
  
  # fishery size comp
  fish_size_lk = 0.0
  offset_fs = 0.0
  for(t in 1:Tfs) {
    # for(l in 1:L) {
    offset_fs = offset_fs - fish_size_iss[t] * sum((fish_size_obs[,t] + g) * log(fish_size_obs[,t] + g))
    fish_size_lk = fish_size_lk - fish_size_iss[t] * sum((fish_size_obs[,t] + g) * log(fish_size_pred[,t] + g))
    # }
  }
  fish_size_lk = fish_size_lk - offset_fs
  
  like_srv = srv_like * srv_wt
  like_fish_age = fish_age_lk * fish_age_wt
  like_srv_age = srv_age_lk * srv_age_wt
  like_fish_size = fish_size_lk * fish_size_wt
  like_rec = sum(log_Rt^2) / (2 * sigmaR^2) + (length(log_Rt) * log(sigmaR))
  f_regularity = wt_fmort_reg * sum(log_Ft^2)
  
  nll = ssqcatch
  nll = nll + like_srv
  nll = nll + like_fish_age
  nll = nll + like_srv_age
  nll = nll + like_fish_size
  nll = nll + like_rec 
  nll = nll + f_regularity 
  nll = nll - nll_M
  nll = nll - nll_q
  nll = nll + sprpen
  
  # reports -------------------
  RTMB::REPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)   
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)    
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)    
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  
  RTMB::REPORT(spawn_bio_proj)
  RTMB::REPORT(catch_proj)
  RTMB::REPORT(catch_ofl_proj)
  RTMB::REPORT(F40_proj)
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}

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
# )

# run ----
f(pars)
obj <- RTMB::MakeADFun(f, pars, hessian = TRUE,
                       map = map)                 

fit <- nlminb(obj$par, obj$fn, obj$gr)
sd <- RTMB::sdreport(obj)
report <- obj$report(obj$env$last.par.best)
report
obj1 <- RTMB::MakeADFun(f, pars, hessian = TRUE) # ,
# map = list(#og_M = factor(NA),
#            log_q = factor(NA)))             
fit1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
sd1 <- RTMB::sdreport(obj1)
report1 <- obj1$report(obj1$env$last.par.best)
report1
report$catch_pred
report$srv_pred
report$like_srv
report$ssqcatch
report$like_fish_age
report$like_srv_age
report$like_fish_size
report$spawn_bio
report$tot_bio
report$M
report$q
report$Bzero 
report$sigmaR
report$log_mean_R
report$B40


report1$catch_pred
report1$srv_pred
report1$like_srv
report1$ssqcatch
report1$like_fish_age
report1$like_srv_age
report1$like_fish_size
report1$spawn_bio
report1$tot_bio
report1$M
report1$q
report1$Bzero * 0.4
report1$B40

sv = c(126633, 138550, 142971, 142962, 138780, 138752, 138889, 139249, 137494, 131288, 120554, 107594, 96254.4, 88136.6, 79734.1 )
data.frame(years=srv_yrs,
           rep = sv,
           one = report$srv_pred) %>% #,
  # two = report1$srv_pred) %>% 
  tidyr::pivot_longer(-years) %>% 
  ggplot(aes(years, value, color = name)) +
  geom_point() +
  geom_line()

data.frame(year = years,
           f1 = c(0.00769708, 0.0319478, 0.0719912, 0.148687, 0.274143, 0.189075, 0.126403, 0.113601, 0.0853571, 
                  0.0530074, 0.092937, 0.0963922, 0.0712692, 0.0640526, 0.0616836, 0.0529102, 0.013308, 0.0100069, 
                  0.0103002, 0.0110622, 0.0185461, 0.0465545, 0.0410088, 0.0104055, 0.00170176, 0.0020561, 0.00373078, 
                  0.00811458, 0.0106727, 0.0113984, 0.0287564, 0.0479452, 0.0289035, 0.0351226, 0.03329, 0.0198442, 
                  0.0174677, 0.0183077, 0.0328278, 0.0207148, 0.0197433, 0.0209938, 0.032455, 0.0294614, 0.0276125, 
                  0.0301481, 0.0254155, 0.0247846, 0.0246121, 0.0249514, 0.022821, 0.0352467, 0.0359497, 0.0334387, 
                  0.0326463, 0.0300291, 0.0168359, 0.0225628, 0.0276158, 0.0252358, 0.0264193, 0.021821),
           f2 = exp(report$log_mean_F + report$log_Ft) * max(report$slx[,1])) %>% 
  tidyr::pivot_longer(-year) %>% 
  ggplot(aes(year, value, color = name)) +
  geom_point() +
  geom_line()

# Create the data frame
data.frame(
  year = c(1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987, 1988, 1989, 1990, 1991, 1992, 1993, 1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022),
  SpBiom = c(52188.8, 50818.6, 48374.4, 44346.4, 37809.9, 28670.5, 23870.7, 21381.7, 19653.2, 18810, 18769.4, 18156.2, 17605.1, 17586.8, 17827, 18322.2, 19280.1, 21386.8, 24115.3, 27377.4, 30988.2, 34471.3, 36896.1, 39406.6, 43217.8, 47775.1, 52650.9, 57460.5, 61722.9, 65441.1, 68867.5, 71045.1, 71853.8, 73869.6, 75129, 75964.4, 77104.1, 77739.9, 77751.6, 76278.6, 75549, 74931.9, 74504.2, 73718.5, 73626.5, 73941.2, 74078.8, 74343.8, 74260.5, 73639, 72361.3, 70623.3, 67525.4, 64081.6, 60669.3, 57319.3, 54230.1, 51975.4, 49547.2, 47014.8, 44737.8, 42555.3),
  spawn_bio = report$spawn_bio, 
  Tot_biom = c(119809, 117471, 113042, 105405, 92522.8, 74079.7, 64940.7, 60647.7, 57800.3, 56866.8, 57816.3, 58579.4, 60282.1, 63951, 68541.7, 73480.8, 78863, 87174.4, 96564.4, 106369, 116187, 125181, 131554, 138341, 147389, 158062, 168707, 178686, 187393, 194853, 200820, 202824, 200351, 199780, 196931, 194645, 194933, 195825, 196825, 196100, 197516, 198877, 199591, 197577, 194994, 191586, 186684, 181586, 175830, 169580, 162969, 156532, 148267, 140180, 132734, 125763, 119466, 114994, 110257, 105483, 101489, 97950.1),
  totbio = report$tot_bio) %>% 
  tidytable::summarise(spawn_bio = (SpBiom - spawn_bio) / ((SpBiom + spawn_bio) / 2)  ,
                       tot = (Tot_biom - totbio) / ((Tot_biom + totbio) / 2), .by = year) %>% 
  tidyr::pivot_longer(-year) %>% 
  ggplot(aes(year, value, color = name)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0, lty = 3) +
  scale_y_continuous(labels = scales::percent_format())

fish_slx = read.table('clipboard')
data.frame(  
  age = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51),
  fish_slx = fish_slx,
  srv_slx = c(0.00787849, 0.0154587, 0.0301107, 0.0578345, 0.108236, 0.193537, 0.321806, 0.484061, 0.649747, 0.785773, 0.878823, 0.93481, 0.965932, 0.982475, 0.991059, 0.995458, 0.997698, 0.998834, 0.99941, 0.999702, 0.999849, 0.999924, 0.999961, 0.99998, 0.99999, 0.999995, 0.999997, 0.999999, 0.999999, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  fslx = report$slx[,1],
  sslx = report$slx[,2]
) %>% 
  mutate(fish = fish_slx-fslx,
         srv = (srv_slx - sslx)) %>% 
  tidyr::pivot_longer(-age) %>% 
  ggplot(aes(age, value, color = name)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0, lty = 3) +
  scale_y_continuous(labels = scales::percent_format())

df$SpBiom

plot(years, report$catch_pred, pch = 19)
lines(years, report1$catch_pred, pch = 19)

plot(srv_yrs, report$srv_pred, pch = 19, ylim = c(0, 170000))
lines(srv_yrs, report1$srv_pred,)

plot(report$catch_pred, report1$catch_pred)
plot(report$srv_pred, report1$srv_pred)
plot(report$spawn_bio, report1$spawn_bio)
plot(report$tot_bio, report1$tot_bio)
cv = sd_M/1 

(log_M-log(mean_M))^2 / (2.* 0.15^2);

hist(rnorm(1000, mean=0.07, sd=0.015) )
report1$M



diagnostic=TRUE, eps=0.1)
fit <- tmbstan(obj1, chains=1)

cores <- parallel::detectCores()-1
options(mc.cores = cores)
fit <- tmbstan(obj, chains=cores, open_progress=FALSE)

## Can also get ESS and Rhat from rstan::monitor

pairs(fit, pars=names(obj1$par))
class(fit)
names(fit)
launch_shinystan(fit)
traceplot(fit, pars=names(obj$par), inc_warmup=TRUE)


dat <- obj1$env$data
length(dat$ages)

names(obj1$env)
names(obj1$env)

