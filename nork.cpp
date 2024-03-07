#include<TMB.hpp>
#include "functions.hpp" // functions for the model - housed them seperate a la an R package
  using namespace density; 


// testing the setup of the AFSC rockfish model in a TMB framework

template<class Type>
Type objective_function<Type>::operator() ()
{
  // inputs --------------------------------
  DATA_INTEGER(rec_age);         // recruitment age - different from 'usual' since rockfish age error extends beyond plus group
  DATA_VECTOR(years);           // years, dim = T
  DATA_VECTOR(length_bins);     // lengths
  DATA_SCALAR(spawn_mo);
  DATA_VECTOR(waa);    // weight at age
  DATA_VECTOR(wt_mature); // weight at age for mature fish
  // catch
  DATA_VECTOR(catch_obs);
  DATA_INTEGER(yr_catchwt);
  DATA_SCALAR(wt0);
  DATA_SCALAR(wt1);
  // survey biomass
  DATA_IVECTOR(srv_ind);   // survey years
  DATA_VECTOR(srv_obs);    // survey observations
  DATA_VECTOR(srv_sd);     // survey observation standard deviation
  DATA_SCALAR(srv_wt);
  // fishery age comps
  DATA_IVECTOR(fish_age_ind);
  DATA_MATRIX(fish_age_obs);
  DATA_VECTOR(fish_age_iss);
  DATA_SCALAR(fish_age_wt);
  // survey age comps
  DATA_IVECTOR(srv_age_ind);
  DATA_MATRIX(srv_age_obs);
  DATA_VECTOR(srv_age_iss);
  DATA_SCALAR(srv_age_wt);
  //fishery length comps
  DATA_IVECTOR(fish_size_ind);
  DATA_MATRIX(fish_size_obs);
  DATA_VECTOR(fish_size_iss);
  DATA_SCALAR(fish_size_wt);
  // transition matrix
  DATA_MATRIX(age_error);
  DATA_MATRIX(size_age);

  // DATA_SCALAR(sigmaR);

  // parameters -----------------------------
  PARAMETER_VECTOR(log_a50);    // age at 50% selectivity for commercial fishery & survey
  PARAMETER_VECTOR(delta);      // spread of selectivity for commercial fishery & survey
  PARAMETER(log_mean_F);        // log mean F
  PARAMETER_VECTOR(log_Ft);     // log F at time t dim(T)
  PARAMETER(log_M);             // log M
  PARAMETER(log_mean_R);        // log mean recruitment
  PARAMETER_VECTOR(log_Rt);     // log recruitment at time t dim(T)
  PARAMETER_VECTOR(init_logRt); // log recruitment at time t dim(A-2)
  PARAMETER(log_q);             // survey catchability

// setup -----------------------------------------------------------
    int A = age_error.rows();   // total ages - modeled in plus group
    int T = years.size();
    vector<Type> a50 = exp(log_a50);
    Type M = exp(log_M);
    int G = log_a50.size();     // number of fleets
      int indC = 0;             // fishing fleet slx index
    Type spawn_adj = pow(exp(-M), ((spawn_mo - 1) / 12));   // adjustment for spawning month
    Type q = exp(log_q);
    // Type sigmaR2 = pow(sigmaR, 2);
    Type nll = 0.0;             // negative log likelihood

// containers -----------------------------------------------------
    matrix<Type> slx(A,G);          // selectivity for all fleets
    slx.setZero();
    matrix<Type> Fat(A,T);          // F for all ages and years
    Fat.setZero();
    matrix<Type> Zat(A,T);          // Z for all ages and years
    Zat.setZero();
    matrix<Type> Sat(A,T);          // S for all ages and years
    Sat.setZero();
    vector<Type> Ft(T);             // F for all years
    Ft.setZero();
    vector<Type> initNat(A);        // initial numbers at age
    initNat.setZero();
    matrix<Type> Nat(A,T);          // numbers at age
    Nat.setZero();
    vector<Type> totbio(T);         // total biomass
    totbio.setZero();
    vector<Type> ssb(T);            // spawning stock biomass
    ssb.setZero();
    Type Bzero = 0.0;               // unfished biomass
    matrix<Type> Cat(A,T);          // catch at age
    Cat.setZero();
    vector<Type> catch_pred(T);     // predicted catch
    catch_pred.setZero();
    vector<Type> pred0(T);          // first portion of predicted catch
    pred0.setZero();
    vector<Type> pred1(T);          // 2nd portion of predicted catch     
    pred1.setZero();  
    vector<Type> srv_pred = srv_obs.size(); // predicted survey biomass
    srv_pred.setZero();
    matrix<Type> srv_age_pred(srv_age_obs.rows(), srv_age_obs.cols());  // predicted survey age composition
    srv_age_pred.setZero();
    matrix<Type> fish_age_pred(fish_age_obs.rows(), fish_age_obs.cols());   // predicted fishery age composition
    fish_age_pred.setZero();
    matrix<Type> fish_size_pred(fish_size_obs.rows(), fish_size_obs.cols());  // predicted fishery size composition
    // Type rec_pen = 0.0; // recruitment penalty
    // Type Ft_pen = 0.0;  // F-dev penalty
// run functions ----------------------------------------------------
    get_selectivity(slx, a50, delta, rec_age);
    get_survival(Fat, Zat, Sat, Ft, log_mean_F, log_Ft, slx, indC, M);
    get_population(initNat, Nat, Bzero, Sat, wt_mature, spawn_adj, log_mean_R, M, log_Rt, init_logRt);
    get_catch(Cat, catch_pred, pred0, pred1, yr_catchwt, years, Nat, Sat, Fat, Zat, waa);
    get_survey(srv_obs, srv_ind, Nat, slx, waa, q, srv_pred);
    get_fish_age(fish_age_ind, Cat, age_error, fish_age_pred);
    get_srv_age(srv_age_ind, Nat, age_error, slx, srv_age_pred); 
    get_fish_size(fish_size_ind, Cat, size_age, fish_size_pred);
    get_biomass(Nat, waa, wt_mature, totbio, ssb);

// likelihoods -------------------------------------------
    Type catch_likelihood = catch_like(wt0, wt1, yr_catchwt, years, catch_obs, pred0, pred1);
    Type srv_likelihood = survey_like(srv_wt, srv_obs, srv_pred, srv_sd);
    Type age_likelihood = comp_like(fish_age_wt, fish_age_obs, fish_age_iss, fish_age_pred);
    Type srv_age_likelihood = comp_like(srv_age_wt, srv_age_obs, srv_age_iss, srv_age_pred);
    Type length_likelihood = comp_like(fish_size_wt, fish_size_obs, fish_size_iss, fish_size_pred);
    nll = catch_likelihood + srv_likelihood + age_likelihood + srv_age_likelihood + length_likelihood;

// penalties ---------------------------------------------
  // Recruitment penalty
  // for(int t=0; t<init_logRt.size(); ++t){
  //   rec_pen += square(init_logRt(t) + sigmaR2 / 2.) / (2.* sigmaR2);
  //   }
  // for(int t=0; t<log_Rt.size(); ++t){
  //   rec_pen += square(log_Rt(t) + sigmaR2 / 2.) / (2.* sigmaR2);
  //   }
  //   rec_pen += (init_logRt.size() + log_Rt.size()) * log(sigmaR);

  // F-dev penalty
  // Ft_pen = square(log_Ft).sum();

  // nll = nll + rec_pen + Ft_pen;

// report ----------------------------------------------
    REPORT(Bzero);
    REPORT(totbio);
    REPORT(ssb);
    REPORT(slx);
    REPORT(Fat);
    REPORT(Zat);
    REPORT(Sat);
    REPORT(Ft);
    REPORT(M);
    REPORT(Nat);
    REPORT(catch_pred);
    REPORT(srv_pred);
    REPORT(q);
    REPORT(a50);
    REPORT(delta);
    REPORT(fish_age_pred);
    REPORT(srv_age_pred);
    REPORT(fish_size_pred);
    REPORT(init_logRt);
    REPORT(log_Rt);
    REPORT(catch_likelihood);
    REPORT(srv_likelihood);
    REPORT(age_likelihood);
    REPORT(srv_age_likelihood);
    REPORT(length_likelihood);

    return nll;

}

