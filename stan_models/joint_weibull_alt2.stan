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
int N_Icens;
vector<lower=0>[N_obs] lifetime_obs;
vector<lower=0>[N_Rcens] lifetime_Rcens;
vector<lower=0>[N_Icens] lifetime_Icens_Upper;
vector<lower=0>[N_Icens] lifetime_Icens_Lower;
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
vector<lower = lifetime_Rcens>[N_Rcens] y_rcens;
vector<lower = lifetime_Icens_Lower, upper = lifetime_Icens_Upper>[N_Icens] y_ltrunc;
}
transformed parameters {
real<lower = 0> beta;
real<lower = 0> eta;
vector<lower=0>[N_Icens] t_ltrunc;

// calculate Weibull paramaters based on the
// draws from the CDF at t1 and t2.
beta = (fn(t2CDF) - fn(t1CDF)) / log(t_2 / t_1);
eta = exp(log(t_1) - (fn(t1CDF) / beta));

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

// Prior models
// t1CDF ~ beta(par1_t1, par2_t1);
// t2CDF ~ beta(par1_t2, par2_t2);
t1CDF ~ normal(t1_mean, t1_var);
t2CDF ~ normal(t2_mean, t2_var);
}
