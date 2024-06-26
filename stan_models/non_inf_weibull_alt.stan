data {
int N_obs;
int N_Icens;
int N_Rcens;
vector<lower=0>[N_obs] lifetime_obs;
vector<lower=0>[N_Rcens] lifetime_Rcens;
vector<lower=0>[N_Icens] lifetime_Icens_Upper;
vector<lower=0>[N_Icens] lifetime_Icens_Lower;
real<lower = 0> eta_mean;
real<lower = 0> eta_sd;
real<lower = 0> beta_mean;
real<lower = 0> beta_sd;
}
parameters {
real<lower = 0> beta;
real<lower = 0> eta;
vector<lower = lifetime_Rcens>[N_Rcens] Y_Rcens;
vector<lower=lifetime_Icens_Lower, upper=lifetime_Icens_Upper>[N_Icens] Y_Icens;
}
model{
// Data model
// non-censored portion
lifetime_obs ~ weibull(beta, eta);
// right censored portion
Y_Rcens ~ weibull(beta, eta);
// interval portion
Y_Icens ~ weibull(beta, eta);
  
// Prior model
eta ~ normal(eta_mean, eta_sd);
beta ~ normal(beta_mean, beta_sd);
}
