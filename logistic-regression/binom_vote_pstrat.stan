data {
  int<lower=1> N;
  int<lower=1> J_age;
  int<lower=1, upper=J_age> idx_age[N];
  int<lower=1> J_state;
  int<lower=1, upper=J_state> idx_state[N];
  int<lower=1> K[N];
  int<lower=0, upper=max(K)> y[N];
  
  // pstrat params
  int<lower=1> N_ps;
  int<lower=1> cell_elig[N_ps];
  int<lower=1> idx_age_ps[N_ps];
  int<lower=-1> idx_state_ps[N_ps];
}
transformed data {
  vector[50] state_elig = rep_vector(0, 50);
  vector[J_age] n_age = rep_vector(0, J_age);
  
  for (n in 1:N_ps)
    state_elig[idx_state_ps[n]] = state_elig[idx_state_ps[n]] + cell_elig[n];
  for (n in 1:N)
    n_age[idx_age[n]] = n_age[idx_age[n]] + K[n];
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
  int y_rep[N];
  vector[N] p_rep;
  vector[50] N_state_votes = rep_vector(0, 50);
  vector[50] state_prop;
  vector[J_age] age_prop = rep_vector(0, J_age);
  for (n in 1:N) {
    real logit_mu = alpha_age[idx_age[n]] + alpha_state[idx_state[n]];
    y_rep[n] = binomial_rng(K[n], inv_logit(logit_mu));
    p_rep[n] = 1.0 * y_rep[n] / K[n];
    age_prop[idx_age[n]] = age_prop[idx_age[n]] + y_rep[n];
  }
  age_prop = age_prop ./ n_age;
  for (n in 1:N_ps) {
    int idx_state_n = idx_state_ps[n];
    real a_state = idx_state_n > J_state ? normal_rng(0, sigma_state) : alpha_state[idx_state_n];
    real ps_mu = inv_logit(alpha_age[idx_age_ps[n]] + a_state) * cell_elig[n];
    N_state_votes[idx_state_n] = N_state_votes[idx_state_n] + ps_mu;
  }
  state_prop = N_state_votes ./ state_elig;
}
