data {
int N_obs;
int N_Rcens;
int N_Icens;
real lifetime_obs[N_obs];
real lifetime_Rcens[N_Rcens];
real lifetime_Icens_Upper[N_Icens];
real lifetime_Icens_Lower[N_Icens];
}
parameters {
real<lower = 0> beta;
real<lower = 0> eta;
}
model{
// Data model
// non-censored portion
for(i in 1:N_obs){
  target += weibull_lpdf(lifetime_obs[i]|beta, eta);
}
// right censored portion
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

// Prior model
eta ~ normal(1000, 500);
beta ~ normal(1.1, 0.5);
}
