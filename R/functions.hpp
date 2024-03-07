/*
  t,...,T - time steps
  a,...,A - age classes
  l,...,L - length bins
  g,...,G - fleets
  
*/

template <class Type>
Type square(Type x){return x*x;}

// base selectivity function
template<class Type>
Type sel(int age, Type a50, Type delta, int first_age) {
     return 1./(1. + exp(-log(19) * (age + first_age - a50) / delta));
  }

// run selectivity function for all ages and fleets
template <class Type>
void get_selectivity(matrix<Type>& slx, vector<Type>& a50, vector<Type>& delta, int first_age) {
  int A = slx.rows(); // number of age classes
  int G = slx.cols(); // number of fleets
  for(int a = 0; a < A; a++){
      for(int g = 0; g < G; g++){
        slx(a,g) = sel(a, a50(g), delta(g), first_age); // selectivity for each fleet
      }
  }
}

// get mortality rates
template<class Type>
void get_survival(matrix<Type>& Fat, matrix<Type>& Zat, matrix<Type>& Sat, vector<Type>& Ft, 
                    Type& log_mean_F, vector<Type>& log_Ft, matrix<Type>& slx, int indC, Type& M) {
    int T = Ft.size();

    for(int t = 0; t < T; ++t) {
        Ft(t) = exp(log_mean_F + log_Ft(t)); // Fishing mortality
        Fat.col(t) = Ft(t) * slx.col(indC);  // Fully selected fishing mortality
        Zat.col(t) = (Fat.col(t).matrix().array() + M).matrix();         // Total mortality  
        Sat.col(t) = exp(-Zat.col(t).array());       // Survivorship
    }
}

// numbers at age 
template<class Type>
void get_population(vector<Type>& initNat, matrix<Type>& Nat, Type& Bzero, const matrix<Type>& Sat, 
                    const vector<Type>& wt_mature, Type spawn_adj, Type& log_mean_R, Type& M, 
                    vector<Type>& log_Rt, vector<Type>& init_logRt) {
    int T = Nat.cols();     // number of time steps
    int A = Nat.rows();     // number of age classes

    initNat(0) = exp(log_mean_R);
    for (int a=1; a<A-1; a++) {
        initNat(a) = (exp(log_mean_R - M * (a-1) + init_logRt(a-1)));   // initial numbers at age
    }
    initNat(A-1) = initNat(A-2) * exp(-M) / (1.0 - exp(-M)) ;

    Bzero = (initNat.array() * wt_mature.array()).sum() * spawn_adj;    // unfished biomass
    
    Nat.col(0) = initNat;
    Nat.row(0) = exp(log_mean_R + log_Rt);
    for (int t = 1; t < T; t++) {
          for(int a=1; a<A-1; a++) {
            Nat(a,t) = Nat(a-1,t-1) * Sat(a-1,t-1);
          }
            Nat(A-1,t) = Nat(A-2,t-1) * Sat(A-2,t-1) + Nat(A-1,t-1) * Sat(A-1,t-1); 
        }
}

// numbers at age 
template<class Type>
void get_population2(vector<Type>& initNat, matrix<Type>& Nat, Type& Bzero, const matrix<Type>& Sat, 
                    const vector<Type>& wt_mature, Type spawn_adj, Type& log_mean_R, Type& M, 
                    vector<Type>& log_Rt, vector<Type>& init_logRt) {
    int T = Nat.cols();     // number of time steps
    int A = Nat.rows();     // number of age classes

    initNat(0) = exp(log_mean_R);
    for (int a=1; a<A-1; a++) {
        initNat(a) = (exp(log_mean_R - M * (a-1) + init_logRt(a-1)));   // initial numbers at age
    }
    initNat(A-1) = initNat(A-2) * exp(-M) / (1.0 - exp(-M)) ;

    Bzero = (initNat.array() * wt_mature.array()).sum() * spawn_adj;    // unfished biomass
    
    Nat.col(0) = initNat;
    Nat.row(0) = exp(log_mean_R + log_Rt);
    for (int t = 1; t < T; t++) {
          for(int a=1; a<A-1; a++) {
            Nat(a,t) = Nat(a-1,t-1) * Sat(a-1,t-1);
          }
            Nat(A-1,t) = Nat(A-2,t-1) * Sat(A-2,t-1) + Nat(A-1,t-1) * Sat(A-1,t-1); 
        }
}

template<class Type>
void get_catch(matrix<Type>& Cat, vector<Type>& catch_pred, vector<Type>& pred0, vector<Type>& pred1, 
                int yr_catchwt, vector<Type>& years, matrix<Type> Nat, matrix<Type> Sat, 
                matrix<Type> Fat, matrix<Type> Zat, vector<Type>& waa){

    int T = Nat.cols();

    vector<Type> catch_ind(T);

    for(int t=0; t<T; t++) {
      catch_ind(t) = years(t) <= yr_catchwt ? 0 : 1;
    }
  
    pred0.setZero();  // set to zero
    pred1.setZero();
    catch_pred.setZero();

    Cat = (Fat.array() / Zat.array() * Nat.array() * (1.0 - Sat.array())).matrix(); // catch at age

    for(int t=0; t<T; t++){
        pred0(t) = (catch_ind(t) == 0) ? (Cat.col(t).array() * waa.array()).sum() : 0;
        pred1(t) = (catch_ind(t) == 1) ? (Cat.col(t).array() * waa.array()).sum() : 0;
        catch_pred(t) = pred0(t) + pred1(t);
    }
}

// catch likelihood for the 2 time blocks
template<class Type>
Type catch_like(Type wt0, Type wt1, int yr_catchwt, vector<Type> years, vector<Type> catch_obs, vector<Type>& pred0, vector<Type>& pred1){
    Type ssqcatch = 0.0;
    Type g = 0.00001;
    int T = catch_obs.size();
    vector<Type> catch_ind(T);

    for(int t=0; t<T; t++) {
      catch_ind(t) = years(t) <= yr_catchwt ? 0 : 1;
    }
    
    vector<Type> catch_obs0 = (catch_ind.array() == 0).select(catch_obs, 0);
    vector<Type> catch_obs1 = (catch_ind.array() == 1).select(catch_obs, 0);

    ssqcatch += wt0 * (log(catch_obs0.array() + g) - log(pred0.array() + g)).square().sum();
    ssqcatch += wt1 * (log(catch_obs1.array() + g) - log(pred1.array() + g)).square().sum();

    return(ssqcatch);
}

template<class Type>
void get_survey(vector<Type> obs, vector<int> ind, matrix<Type> nat, matrix<Type> slx, vector<Type> waa, Type& q, vector<Type>& srv_pred) {
    int T = nat.cols();
    int A = nat.rows();
    vector<Type> pred(T);
    pred.setZero();
    int isrv = 0;  

    for(int t=0; t<T; t++){
      if(ind(t) == 1) {
        for(int a=0; a<A; a++){
            srv_pred(isrv) += nat(a,t) * slx(a,1) * waa(a);
        }
        pred(t) *= q;
        ++isrv;
      }
    }
}

template<class Type>
Type survey_like(Type wt, vector<Type>& obs, vector<Type>& pred, vector<Type>& sd) {
    int T = obs.size();

    Type like = 0.0;
    for(int t=0; t<T; t++){
        like += wt * dnorm(log(obs(t)), log(pred(t)), sd(t), true);
    }
    return(like);
}


template<class Type>
void get_fish_age(vector<int>& ind, matrix<Type>& cat, matrix<Type>& age_error, matrix<Type>& fish_age_pred) {
    int T = cat.cols();
    int icomp = 0;  
    
    for(int t=0; t<T; ++t) {
        if(ind(t) == 1) {
          vector<Type> temp_cat = cat.col(t) / sum(vector<Type>(cat.col(t)));
          vector<Type> temp_ae = (temp_cat.matrix().transpose()) * age_error;
          fish_age_pred.col(icomp) = temp_ae / temp_ae.sum();
          ++icomp;
        }
    }
}

template<class Type>
void get_srv_age(vector<int>& ind, matrix<Type>& nat, matrix<Type>& age_error, matrix<Type>& slx, matrix<Type>& srv_age_pred) {
    int T = nat.cols();
    int icomp = 0;  
    
    for(int t=0; t<T; ++t) {
        if(ind(t) == 1) {
          vector<Type> temp_nat = (nat.col(t).array() * slx.col(1).array()) / (nat.col(t).array() * slx.col(1).array()).sum();
          vector<Type> temp_ae = (temp_nat.matrix().transpose()) * age_error;
          srv_age_pred.col(icomp) = temp_ae / temp_ae.sum();
          ++icomp;
        }
    }
}

template<class Type>
void get_fish_size(vector<int>& ind, matrix<Type>& cat, matrix<Type>& sizeage, matrix<Type>& fish_size_pred) {
    int T = cat.cols();
    int icomp = 0;  
    for(int t=0; t<T; ++t) {
        if(ind(t) == 1) {
          vector<Type> temp_cat = cat.col(t) / sum(vector<Type>(cat.col(t)));
          vector<Type> temp_sa = (temp_cat.matrix().transpose() * sizeage).transpose();
          fish_size_pred.col(icomp) = temp_sa / temp_sa.sum();
          ++icomp;
        }
    }
}

template<class Type>
Type comp_like(Type& wt, matrix<Type>& obs, vector<Type>& iss, matrix<Type>& pred) {
    int T = obs.cols();
    int A = obs.rows();
    Type g = 0.00001;
    Type offset = 0.0;
    Type like = 0.0;
    vector<Type> temp(A);  

    for(int t=0; t<T; t++){
      temp = obs.col(t) / iss(t);
      offset -= iss(t) * ((obs.col(t).array() + g) * log(pred.col(t).array() + g)).sum();
      like -= iss(t) * (vector<Type>((temp.array() + g) * log(pred.col(t).array() + g))).sum();
      like -= wt * offset;
    }
    return(like);
}

template<class Type>
void get_biomass(const matrix<Type>& Nat, const vector<Type>& waa, const vector<Type>& wt_mature,
                        vector<Type>& totbio, vector<Type>& ssb) {
    int T = Nat.cols();
    totbio.resize(T);
    ssb.resize(T);

    for(int t = 0; t < T; ++t) {
        totbio(t) = (Nat.col(t).array() * waa).sum();
        ssb(t) = (Nat.col(t).array() * wt_mature).sum();
    }
}

