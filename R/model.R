f <- function(pars) {
  require(RTMB)
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
  
  # priors -----------------
  # prefer using dnorm but want to match admb exactly
  # won't affect nll location, but values will be slightly different
  # admb priors do not include a constant that is in dnorm
  # nll_M = dnorm(M, mean_M, cv_M, TRUE)
  # nll_q = dnorm(q, mean_q, cv_q, TRUE)
  # nll_sigmaR = dnorm(sigmaR, mean_sigmaR, cv_sigmaR, TRUE)
  
  nll_M = (log_M - log(mean_M))^2 / (2 * cv_M^2)
  nll_q = (log_q - log(mean_q))^2 / (2 * cv_q^2)
  nll_sigmaR= (log(sigmaR / mean_sigmaR))^2 / (2 * cv_sigmaR^2)
  
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
  
  # survey biomass ----
  isrv = 1
  srv_like = 0.0
  # survey biomass & likelihood
  for(t in 1:T) {
    if(srv_ind[t]==1) {
      # for(a in 1:A) {
      srv_pred[isrv] = sum(Nat[,t] * slx[,2] * waa)
      # }
      srv_pred[isrv] = srv_pred[isrv] * q
      # srv_like = srv_like + sum((log(srv_obs[isrv]) - log(srv_pred[isrv]))^2 /
      #                             (2 * (srv_sd[isrv] / srv_obs[isrv])^2))
      srv_like = srv_like + sum((srv_obs[isrv]-srv_pred[isrv])^2/ (2.*(srv_sd[isrv]^2)))
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
  data.frame(log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # filter years 1979:
  df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
  pred_rec = mean(df$pred_rec)
  stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1))
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
  # like_rec = sum(log_Rt^2) / (2 * sigmaR^2) + (length(log_Rt) * log(sigmaR))
  like_rec = sum((c(log_Rt , init_log_Rt) + sigmaR * sigmaR / 2)^2) / (2 * sigmaR^2)
  f_regularity = wt_fmort_reg * sum(log_Ft^2)
  
  nll = 0.0
  nll = ssqcatch
  nll = nll + like_srv
  nll = nll + like_fish_age
  nll = nll + like_srv_age
  nll = nll + like_fish_size
  nll = nll + like_rec * wt_rec_var
  nll = nll + f_regularity
  nll = nll - nll_M
  nll = nll - nll_q
  nll = nll + sprpen
  
  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(stdev_rec)
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}
f1 <- function(pars) {
  require(RTMB)
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
  
  # priors -----------------
  nll_M = nll_q = nll_sigmaR = 0.0
  nll_M = dnorm(M, mean_M, sd_M, TRUE)
  nll_q = dnorm(q, mean_q, sd_q, TRUE)
  nll_a50C = dnorm(a50C, mean_a50C, sd_a50C, TRUE)
  nll_deltaC = dnorm(deltaC, mean_deltaC, sd_deltaC)
  nll_a50S = dnorm(a50S, mean_a50S, sd_a50S, TRUE)
  nll_deltaS = dnorm(deltaS, mean_deltaS, sd_deltaS)
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
  
  # survey biomass ----
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
      # srv_like = srv_like + sum((srv_obs[isrv]-srv_pred[isrv])^2/ (2.*(srv_sd[isrv]^2)))
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
  data.frame(log_Rt = log_Rt,
             pred_rec = Nat[1,],
             year = years) -> df
  # filter years 1979:
  df = df[years>=(1977+ages[1]) & years<=(max(years)-ages[1]),]
  pred_rec = mean(df$pred_rec)
  stdev_rec = sqrt(sum((df$log_Rt - mean(df$log_Rt))^2) / (length(df$log_Rt) - 1))
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
  
  nll = 0.0
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
  nll = nll - nll_a50C 
  nll = nll - nll_deltaC 
  nll = nll - nll_a50S 
  nll = nll - nll_deltaS
  
  # reports -------------------
  RTMB::REPORT(ages)
  RTMB::REPORT(years)
  RTMB::REPORT(M)
  RTMB::REPORT(a50C)
  RTMB::REPORT(deltaC)
  RTMB::REPORT(a50S)
  RTMB::REPORT(deltaS)
  RTMB::REPORT(q)
  RTMB::REPORT(sigmaR)
  RTMB::REPORT(log_mean_R)
  RTMB::REPORT(log_Rt)
  RTMB::REPORT(log_mean_F)
  RTMB::REPORT(log_Ft)
  RTMB::REPORT(waa)
  RTMB::REPORT(wt_mature)
  RTMB::REPORT(yield_ratio)
  RTMB::REPORT(Fat)
  RTMB::REPORT(Zat)
  RTMB::REPORT(Sat)
  RTMB::REPORT(Cat)
  RTMB::REPORT(Nat)
  RTMB::REPORT(slx)
  RTMB::REPORT(Ft)
  RTMB::REPORT(catch_pred)
  RTMB::REPORT(srv_pred)
  
  RTMB::REPORT(fish_age_pred)
  RTMB::REPORT(srv_age_pred)
  RTMB::REPORT(fish_size_pred)
  
  RTMB::REPORT(tot_bio)
  RTMB::REPORT(spawn_bio)
  RTMB::REPORT(spawn_fract)
  RTMB::REPORT(B0)
  RTMB::REPORT(B40)
  RTMB::REPORT(B35)
  RTMB::REPORT(F35)
  RTMB::REPORT(F40)
  RTMB::REPORT(F50)
  RTMB::REPORT(pred_rec)
  RTMB::REPORT(stdev_rec)
  RTMB::REPORT(ssqcatch)
  RTMB::REPORT(like_srv)
  RTMB::REPORT(like_fish_age)
  RTMB::REPORT(like_srv_age)
  RTMB::REPORT(like_fish_size)
  RTMB::REPORT(like_rec)
  RTMB::REPORT(f_regularity)
  RTMB::REPORT(sprpen)
  RTMB::REPORT(nll_q)
  RTMB::REPORT(nll_M)
  RTMB::REPORT(nll)
  # nll = 0.0
  return(nll)
}
