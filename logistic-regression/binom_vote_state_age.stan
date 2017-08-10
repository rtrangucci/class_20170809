data {
  int<lower=1> N;
  int<lower=1> J_age;
  int<lower=1, upper=J_age> idx_age[N];
  int<lower=1> J_state;
  int<lower=1, upper=J_state> idx_state[N];
  int<lower=1> K[N];
  int<lower=0, upper=max(K)> y[N];
}
parameters {
  real alpha;
  vector[J_age] eta_age;
  vector[J_state] eta_state;
  real<lower=0> sigma_age;
  real<lower=0> sigma_state;
}
transformed parameters {
  vector[J_age] alpha_age;
  vector[J_state] alpha_state;
  
  alpha_age = alpha + sigma_age * eta_age; // why this intercept to add alpha?
  alpha_state = sigma_state * eta_state;
}
model {
  // priors
  alpha ~ normal(0.7, 1);
  eta_age ~ normal(0, 1);
  eta_state ~ normal(0, 1);
  
  sigma_age ~ normal(0, 1);
  sigma_state ~ normal(0, 1);
  
  y ~ binomial_logit(K, alpha_age[idx_age] + alpha_state[idx_state]);
}
generated quantities {
}