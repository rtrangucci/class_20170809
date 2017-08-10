data {
  int<lower=1> N;
  int<lower=1> D;
  int<lower=1> J_age;
  int<lower=1, upper=J_age> idx_age[N];
  int<lower=1> J_state;
  int<lower=1, upper=J_state> idx_state[N];
  int<lower=1> J_edu;
  int<lower=1, upper=J_state> idx_edu[N];
  int<lower=1> K[N];
  int<lower=0, upper=max(K)> y[N];
  
  // QR components
  matrix[D,D] R_inv; 
  matrix[N,D] Q_scale;
  row_vector[D] c_means;
  
  // pstrat params
  int<lower=1> N_ps;
  int<lower=1> cell_elig[N_ps];
  int<lower=1> idx_age_ps[N_ps];
  int<lower=1> idx_state_ps[N_ps];
  int<lower=1> idx_edu_ps[N_ps];
  vector[N_ps] age_vec_ps;
  vector[N_ps] state_pres_vote_ps;
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
  vector[J_edu] eta_edu;
  vector[D] gamma;
  real<lower=0> scale;
  simplex[3] prop_var;
}
transformed parameters {
  vector[J_age] alpha_age;
  vector[J_state] alpha_state;
  vector[J_edu] alpha_edu;
  real sigma_age;
  real sigma_state;
  real sigma_edu;
  
  {
    vector[3] sigmas;
    sigmas = 3 * square(scale) * prop_var;
    sigma_age = sigmas[1];
    sigma_state = sigmas[2];
    sigma_edu = sigmas[3];
  }
  alpha_age = sigma_age * eta_age;
  alpha_state = sigma_state * eta_state;
  alpha_edu = sigma_edu * eta_edu;
}
model {
  vector[N] mu;
  // priors
  alpha ~ normal(0.7, 1);
  eta_age ~ normal(0, 1);
  eta_state ~ normal(0, 1);
  eta_edu ~ normal(0, 1);
  gamma ~ normal(0, 1);
  
  scale ~ exponential(1);
  prop_var ~ dirichlet([2, 2, 2]'); // note inline vector creation
  
  y ~ binomial_logit(K, alpha
                        + alpha_age[idx_age]
                        + alpha_state[idx_state]
                        + alpha_edu[idx_edu]
                        + Q_scale * gamma);
}
generated quantities {
  vector[D] beta = R_inv * gamma;
  real alpha_prime = alpha - c_means * beta;
  int y_rep[N];
  vector[N] p_rep;
  vector[50] N_state_votes = rep_vector(0, 50);
  vector[50] state_prop;
  vector[J_age] age_prop = rep_vector(0, J_age);
  for (n in 1:N) {
    real logit_mu = alpha + alpha_age[idx_age[n]] + alpha_state[idx_state[n]]
               + Q_scale[n] * gamma;
    y_rep[n] = binomial_rng(K[n], inv_logit(logit_mu));
    p_rep[n] = 1.0 * y_rep[n] / K[n];
    age_prop[idx_age[n]] = age_prop[idx_age[n]] + y_rep[n];
  }
  age_prop = age_prop ./ n_age;
  for (n in 1:N_ps) {
    int idx_state_n = idx_state_ps[n];
    real a_state = idx_state_n > J_state ? normal_rng(0, sigma_state) : alpha_state[idx_state_n];
    real ps_mu = alpha_prime + alpha_age[idx_age_ps[n]] + a_state
                 + state_pres_vote_ps[n] * beta[1] + age_vec_ps[n] * beta[2] ;
    N_state_votes[idx_state_n] = N_state_votes[idx_state_n] + inv_logit(ps_mu) * cell_elig[n];
  }
  state_prop = N_state_votes ./ state_elig;
}