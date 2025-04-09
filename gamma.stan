data {
  int<lower=1> n; // number of data points
  vector<lower=0>[n] y; // observed data
  real<lower=0> mu; // for prior of alpha
  real<lower=0> sigma; // for prior of alpha
  real<lower=0> a0; // for prior of beta
  real<lower=0> b0; // for prior of beta
}
parameters {
  real<lower=0> a;
  real<lower=0> b;
}
model {
  target += gamma_lupdf(y | a, b);
  target += normal_lupdf(a | mu, sigma);
  target += gamma_lupdf(b | a0, b0);
}
