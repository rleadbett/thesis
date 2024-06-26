functions {
// function to simplify the calculation of eta and beta
real fn(real tCDF) {
  return log(-log1m(tCDF));
}
}
data {
// Define the data
int N_obs;
int N_Rcens;
vector<lower=0>[N_obs] lifetime_obs;
vector<lower=0>[N_Rcens] lifetime_Rcens;
// Define the prior
real t_1;
real t_2;
real t1_mean;
real t1_var;
real t2_mean;
real t2_var;
}
// transformed data {
// real par1_t1;
// real par2_t1;
// real par1_t2;
// real par2_t2;
// 
// Calculate the params of the first beta prior
// par1_t1 = ((t1_mean ^ 2) * (1 - t1_mean) / t1_var) - t1_mean;
// par2_t1 = (par1_t1 / t1_mean) - par1_t1;
// Calculate the params of the second beta prior
// par1_t2 = ((t2_mean ^ 2) * (1 - t2_mean) / t2_var) - t2_mean;
// par2_t2 = (par1_t2 / t2_mean) - par1_t2;
// }
parameters {
real<lower = 0, upper = 1> t1CDF;
real<lower = t1CDF, upper = 1> t2CDF;
vector<lower = lifetime_Rcens>[N_Rcens] Y_Rcens;
}
transformed parameters {
real<lower = 0> beta;
real<lower = 0> eta;

// calculate Weibull paramaters based on the
// draws from the CDF at t1 and t2.
beta = (fn(t2CDF) - fn(t1CDF)) / log(t_2 / t_1);
eta = exp(log(t_1) - (fn(t1CDF) / beta));
}
model{
// Data model
// non-censored portion
lifetime_obs ~ weibull(beta, eta);
// right censored portion
Y_Rcens ~ weibull(beta, eta);

// Prior models
// t1CDF ~ beta(par1_t1, par2_t1);
// t2CDF ~ beta(par1_t2, par2_t2);
t1CDF ~ normal(t1_mean, t1_var);
t2CDF ~ normal(t2_mean, t2_var);
}
