data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1, upper=J> idx_J[N];
  vector[N] y;
  vector[N] X;
}
parameters {
  real alpha;
  vector[J] eta;
  real gamma;
  vector[J] beta_std;
  real<lower=0> sigma_y;
  real<lower=0> sigma_theta;
  real<lower=0> sigma_beta;
}
transformed parameters {
  vector[J] theta;
  vector[J] beta;
  
  theta = sigma_theta * eta;
  beta = gamma + sigma_beta * beta_std;
}
model {
  y ~ normal(alpha + theta[idx_J] + X .* beta[idx_J], sigma_y);
  
  sigma_y ~ cauchy(0, 10);
  sigma_theta ~ normal(0, 3);
  alpha ~ normal(100, 20);
  eta ~ normal(0, 1); 
  beta_std ~ normal(0, 1);
  sigma_beta ~ normal(0, 1);
  gamma ~ normal(0, 10);
}
generated quantities {
  vector[N] y_rep;
  
  for (n in 1:N)
    y_rep[n] = normal_rng(alpha + theta[idx_J[n]] + X[n] * beta[idx_J[n]],
                          sigma_y);
}