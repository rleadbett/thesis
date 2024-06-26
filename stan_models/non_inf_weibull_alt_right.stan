data {
int N_obs;
int N_Rcens;
vector<lower=0>[N_obs] lifetime_obs;
vector<lower=0>[N_Rcens] lifetime_Rcens;
real<lower = 0> eta_mean;
real<lower = 0> eta_sd;
real<lower = 0> beta_mean;
real<lower = 0> beta_sd;
}
parameters {
real<lower = 0> beta;
real<lower = 0> eta;
vector<lower = lifetime_Rcens>[N_Rcens] Y_Rcens;
}
model{
// Data model
// non-censored portion
lifetime_obs ~ weibull(beta, eta);
// right censored portion
Y_Rcens ~ weibull(beta, eta);
  
// Prior model
eta ~ normal(eta_mean, eta_sd);
beta ~ normal(beta_mean, beta_sd);
}
