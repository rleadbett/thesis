data {
int N_obs;
int N_Rcens;
int N_Icens;
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
vector<lower = lifetime_Rcens>[N_Rcens] y_rcens;
vector<lower = lifetime_Icens_Lower, upper = lifetime_Icens_Upper>[N_Icens] y_ltrunc;
}
transformed parameters {
vector<lower=0>[N_Icens] t_ltrunc;

// calculate the imputed truncation time
t_ltrunc = y_ltrunc - lifetime_Icens_Lower;
}
model{
// Data model
// the likelihood model for the imputed data, accounting for left truncation
lifetime_obs ~ weibull(beta, eta);
y_rcens ~ weibull(beta, eta);
for (i in 1:N_Icens) {
  y_ltrunc[i] ~ weibull(beta, eta) T[t_ltrunc[i], ];
}

  
// Prior model
eta ~ normal(eta_mean, eta_sd);
beta ~ normal(beta_mean, beta_sd);
}
