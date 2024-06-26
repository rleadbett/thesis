// Joint-informative-prior
//
// This is a Bayesian model for censored lifetime data 
// which followes a two-paramater Weibull distribution.
// An informative joint prior on the shape and scale 
// paramaters of the weibull distribution (Kaminskiy2017).
//

functions {
  real fn(real tCDF) {
    return log(-log1m(tCDF));
  }
}

data {
  int N_obs;
  int N_Icens;
  int N_Rcens;
  real lifetime_obs[N_obs];
  real lifetime_Rcens[N_Rcens];
  real lifetime_Icens_Upper[N_Icens];
  real lifetime_Icens_Lower[N_Icens];
  real t_1;   // should be 3.82
  real t_2;   // should be 15
  real t1_mean;
  real t1_var;
  real t2_mean;
  real t2_var;
}
transformed data {
  //real par1_t1;
  //real par2_t1;
  //real par1_t2;
  //real par2_t2;

  // Calculate the params of the first beta prior
  //par1_t1 = ((t1_mean ^ 2) * (1 - t1_mean) / t1_var) - t1_mean;
  //par2_t1 = (par1_t1 / t1_mean) - par1_t1;
  // Calculate the params of the second beta prior
  //par1_t2 = ((t2_mean ^ 2) * (1 - t2_mean) / t2_var) - t2_mean;
  //par2_t2 = (par1_t2 / t2_mean) - par1_t2;
}
parameters {
  real<lower = 0, upper = 1> t1CDF;
  real<lower = t1CDF, upper = 1> t2CDF;  // the CDF at t2 must be greater than at t1
}

transformed parameters {
  real<lower = 0> beta;
  real<lower = 0> eta;

  // calculate Weibull paramaters based on the
  // draws from the CDF at t1 and t2.
  beta = (fn(t2CDF) - fn(t1CDF)) / log(t_2 / t_1);
  eta = exp(log(t_1) - (fn(t1CDF) / beta));
}

model {
  // Likelihood
  // non-censored portion
  for(i in 1:N_obs){
    target += weibull_lpdf(lifetime_obs[i]|beta, eta);
  }
  // censored portion
  for(j in 1:N_Rcens){
    target += weibull_lccdf(lifetime_Rcens[j]|beta, eta);
  }
  // interval portion
  for(k in 1:N_Icens){
    target += log_diff_exp(
      weibull_lcdf(lifetime_Icens_Upper[k]|beta, eta),
      weibull_lcdf(lifetime_Icens_Lower[k]|beta, eta)
    );
  }
  
  // Prior models
  // The prior was constructed by simulateing 100 datasets of size 
  // n = 100 from the true Weibull distribution and estimating the 
  // paramaters via MLE and calculating to value of the estimated 
  // CDF at t1 and t2 to get a distribution.
  t1CDF ~ normal(t1_mean, t1_var);//beta(par1_t1, par2_t1);
  t2CDF ~ normal(t2_mean, t2_var);//beta(par1_t2, par2_t2);
}
