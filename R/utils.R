get_survival <- function(post, obj, reps) {
  
  mapit = function(i, post, obj) {
    surv_i = obj$report(post[i,-ncol(post)])$Nat[,T] * obj$report(post[i,-ncol(post)])$Sat[,T]
    list(surv_i = surv_i)
  }
  results = purrr::map(1:reps, mapit, post=post, obj=obj)
  do.call(cbind, map(results, "surv_i"))
}

# obj$report(post[i,-ncol(post)])$q

get_bio <- function(post, obj, reps) {
  
  mapit = function(i, post, obj) {
    surv_i = obj$report(post[i,-ncol(post)])$spawn_bio
    list(surv_i = surv_i)
  }
  results = purrr::map(1:reps, mapit, post=post, obj=obj)
  do.call(cbind, map(results, "surv_i")) 
    
}
# post = post[200:1000,]
# bb = get_bio(post, obj1, nrow(post))
# bb %>% 
#   as.data.frame() %>% 
#     mutate(year = years) %>% 
#   pivot_longer(-year) %>% 
#   ggplot(aes(year, value, group = name)) + 
#   geom_line(alpha = 0.2) +
#   geom_hline(yintercept = c(report1$B40, report1$B35))
  

proj_bio <- function(report, obj=NULL, post=NULL, reps = 500) {
  
  # values
  F40 = report$F40
  F35 = report$F35
  B40 = report$B40
  Nat = report$Nat
  Sat = report$Sat
  slx = report$slx
  ages = report$ages
  years = report$years
  waa = report$waa
  wt_mature = report$wt_mature
  spawn_frac = report$spawn_fract
  yield_ratio = report$yield_ratio
  M = report$M
  pred_rec = report$pred_rec
  stdev_rec = report$stdev_rec 
  A = nrow(Nat) # number of ages
  T = ncol(Nat) # number of years
  
  # storage
  if(!is.null(post)) {
    surv = get_survival(post=post, obj=obj1, reps = reps)
    Tproj = 15
    N = Fabc_tot = Fofl_tot = replicate(reps, matrix(0, A, Tproj), simplify = FALSE)
    F40_proj = F35_proj = spawn_bio = tot_bio = replicate(reps, rep(0, Tproj), simplify = FALSE)
    for(i in 1:reps) {
      F40_proj[[i]] = rep(F40, Tproj)
      F35_proj[[i]] = rep(F35, Tproj)
      
      # total F
      Fabc_tot[[i]][,1] = report$slx[,1] * F40
      Fofl_tot[[i]][,1] = report$slx[,1] * F35
      
      # populate abundance
      N[[i]][1,] = exp(log(pred_rec) - (stdev_rec)^2 / 2 + stdev_rec + rnorm(15))
      
      for(a in 1:(A-1)) {
        N[[i]][a+1,1] = surv[,i][a]
      }
      N[[i]][A,1] = surv[,i][A-1] + surv[,i][A]
      spawn_bio[i][1] = sum(N[[i]][,1] * exp(-yield_ratio * unlist(Fabc_tot[[i]][,1]) - M)^spawn_frac * wt_mature)
      tot_bio[[i]][1] = sum(N[[i]][,1] * waa)   
      
      for(t in 1:Tproj) {
        # tier check
        if((spawn_bio[[i]][t] / B40) <= 1) {
          F40_proj[[i]][t] = F40_proj[[i]][t] * (spawn_bio[[i]][t] / B40 - 0.05) / 0.95
          F35_proj[[i]][t] = F35_proj[[i]][t] * (spawn_bio[[i]][t] / B40 - 0.05) / 0.95
        }
        # update 
        Fabc_tot[[i]][,t] = report$slx[,1] * F40_proj[[i]][t]
        Fofl_tot[[i]][,t] = report$slx[,1] * F35_proj[[i]][t]
        Z_proj = unlist(Fabc_tot) + M
        Zofl_proj = unlist(Fofl_tot) + M
        S_proj = exp(-Z_proj)
        
        # catch
        # Cat_proj[,t] = yield_ratio * N_proj[,t] * Fabc_tot_proj / Z_proj * (1 - S_proj)
        # Cat_ofl_proj[,t] = yield_ratio * N_proj[,t] * Fofl_tot_proj / Zofl_proj * (1 - exp(-Zofl_proj))
        
        if(t<Tproj) {
          for(a in 1:(A-1)){
            N[[i]][a+1,t+1] = N[[i]][a,t] * exp(-yield_ratio * unlist(Fabc_tot[i])[a] - M)
          }
          N[[i]][A,t+1] = N[[i]][A-1,t] * exp(-yield_ratio * unlist(Fabc_tot[i])[A-1] - M) +
            N[[i]][A,t] * exp(-yield_ratio * unlist(Fabc_tot[i])[A] - M)
          
          tot_bio[[i]][t+1] = sum(N[[i]][,t+1] * waa)  
          spawn_bio[[i]][t+1] = sum(N[[i]][,t+1] * exp(-yield_ratio * unlist(Fabc_tot[[i]][,t]) - M)^spawn_frac * wt_mature)
        }
      }
    }
    
    do.call(rbind, spawn_bio) %>% 
      as.data.frame() %>% 
      mutate(sim = 1:n()) %>% 
      tidyr::pivot_longer(-sim) %>% 
      mutate(year = as.numeric(gsub('V','', name)) + max(years),
             id = 'spawn_bio') %>% 
      select(-name) %>% 
      dplyr::bind_rows(
        do.call(rbind, tot_bio) %>% 
          as.data.frame() %>% 
          mutate(sim = 1:n()) %>% 
          tidyr::pivot_longer(-sim) %>% 
          mutate(year = as.numeric(gsub('V','', name)) + max(years),
                 id = 'tot_bio') %>% 
          select(-name))
    
  } else {
    Tproj = 2
    N = Cat = Cat_ofl= Zabc = Zofl = S = matrix(0, A, Tproj)
    tot_bio = spawn_bio = F40_proj = F35_proj= rep(0, Tproj)
    # setup
    F40_proj[1] = F40
    F35_proj[1] = F35
    
    # total F
    Fabc_tot = slx[,1] * F40_proj[1]
    Fofl_tot = slx[,1] * F35_proj[1]
    
    # first projection year
    N[1,] = pred_rec
    for(a in 1:(A-1)) {
      N[a+1,1] = Nat[a,T] * Sat[a,T]
    }
    N[A,1] = Nat[A-1,T] * Sat[A-1,T] + Nat[A,T] * Sat[A,T]
    spawn_bio[1] = sum(N[,1] * exp(-yield_ratio * Fabc_tot - M)^spawn_frac * wt_mature)
    
    for(t in 1:Tproj) {
      # tier check
      if((spawn_bio[t] / B40) > 1) {
        F40_proj[t] = F40
        F35_proj[t] = F35
      } else {
        F40_proj[t] = F40_proj[t] * (spawn_bio[t] / B40 - 0.05) / 0.95
        F35_proj[t] = F35_proj[t] * (spawn_bio[t] / B40 - 0.05) / 0.95
      }
      # update
      Fabc_tot = report$slx[,1] * F40_proj[t]
      Fofl_tot = report$slx[,1] * F35_proj[t]
      Z = Fabc_tot + M
      Zofl = Fofl_tot + M
      S = exp(-Z)
      
      # catch
      Cat[,t] = yield_ratio * N[,t] * Fabc_tot / Z * (1 - S)
      Cat_ofl[,t] = yield_ratio * N[,t] * Fofl_tot / Zofl * (1 - exp(-Zofl))
      
      if(t<Tproj) {
        for(a in 1:(A-1)){
          N[a+1,t+1] = N[a,t] * exp(-yield_ratio * Fabc_tot[a] - M)
        }
        N[A,t+1] = N[A-1,t] * exp(-yield_ratio * Fabc_tot[A-1] - M) +
          N[A,t] * exp(-yield_ratio * Fabc_tot[A] - M)
        
        tot_bio[t+1] = sum(N[,t+1] * waa)
        spawn_bio[t+1] = sum(N[,t+1] * exp(-yield_ratio * Fabc_tot - M)^spawn_frac * wt_mature)
      }
    }
    catch = colSums(Cat * waa / yield_ratio)
    catch_ofl = colSums(Cat_ofl * waa / yield_ratio)
    tot_bio = colSums(N * waa)
    
    data.frame(year = max(years)+1:Tproj,
               spawn_bio = spawn_bio,
               tot_bio = tot_bio,
               catch_abc = catch,
               catch_ofl = catch_ofl,
               F40 = F40_proj,
               F35 = F35_proj)
  }
}

single_likelihood <- function(param_name, param_values, obj, fit) {
  n <- length(param_values)
  log_likelihoods <- numeric(n)
  
  for (i in 1:n) {
    params <- fit$par
    params[param_name] <- param_values[i]
    log_likelihoods[i] <- obj$fn(params)
  }
  
  return(log_likelihoods)
}


# Function to profile the likelihood for multiple parameters
profile_likelihood <- function(param_names, param_ranges, obj, fit) {
  grid <- expand.grid(param_ranges)
  log_likelihoods <- numeric(nrow(grid))
  
  for (i in 1:nrow(grid)) {
    params <- fit$par
    for (j in 1:length(param_names)) {
      params[param_names[j]] <- grid[i, j]
    }
    log_likelihoods[i] <- obj$fn(params)
  }
  
  return(data.frame(grid, log_likelihood=log_likelihoods))
}
