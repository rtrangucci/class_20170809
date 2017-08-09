data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1, upper=J> idx_J[N];
  vector[N] y;
}
parameters {
  real alpha;
  vector[J] eta;
  real<lower=0> sigma_y;
  real<lower=0> sigma_theta;
}
transformed parameters {
  vector[J] theta;
  
  theta = sigma_theta * eta;
}
model {
  y ~ normal(alpha + theta[idx_J], sigma_y);
  
  sigma_y ~ cauchy(0, 10);
  sigma_theta ~ normal(0, 3);
  alpha ~ normal(100, 20);
  eta ~ normal(0, 1); 
}