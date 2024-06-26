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
vector<lower = lifetime_Icens_Lower, upper = lifetime_Icens_Upper>[N_Icens] Y_Ltrunc;
}
transformed parameters {
vector<lower = 0>[N_Icens] t_Ltrunc;
// calculate the imputed truncation time
t_Ltrunc = Y_Ltrunc - lifetime_Icens_Lower;
}
model{
// Data model
// the likelihood model for the imputed data, accounting for left truncation
target += weibull_lpdf(lifetime_obs| beta, eta);
target += weibull_lpdf(Y_Rcens| beta, eta);
target += weibull_lpdf(Y_Ltrunc| beta, eta) - weibull_lccdf(t_Ltrunc| beta, eta);

// Impute the partialy observed lifetimes (censored lifetimes)
// right censored portion
Y_Rcens ~ weibull(beta, eta);
// interval portion
Y_Ltrunc ~ weibull(beta, eta);
  
// Prior model
eta ~ normal(eta_mean, eta_sd);
beta ~ normal(beta_mean, beta_sd);
}
