data { 
  int<lower=2> K;
  int<lower=0> N;
  int<lower=1> D;
  int<lower=1,upper=K> y[N];
  matrix[N, D] X;
  int<lower=1> J_age;
  int<lower=1> J_eth;
  int<lower=1> J_edu;
  int<lower=1> G;
  int<lower=1, upper=J_eth> idx_eth[N];
  int<lower=1, upper=J_age> idx_age[N];
  int<lower=1, upper=J_edu> idx_edu[N];
}
transformed data {
  row_vector[D] zeros = rep_row_vector(0, D);
  row_vector[J_age] zeros_age = rep_row_vector(0, J_age);
  row_vector[J_eth] zeros_eth = rep_row_vector(0, J_eth);
  row_vector[J_edu] zeros_edu = rep_row_vector(0, J_edu);
  matrix[D, N] t_X = X';
}
parameters {
  matrix[K - 1, D] beta_raw;
  vector[K - 1] alpha_raw;
  matrix[K - 1, J_age] eta_age;
  matrix[K - 1, J_eth] eta_eth;
  matrix[K - 1, J_edu] eta_edu;
  vector<lower=0>[K - 1] sigma_age;
  vector<lower=0>[K - 1] sigma_eth;
  vector<lower=0>[K - 1] sigma_edu;
  vector<lower=0>[G] sigma_inter_eqn;
}
transformed parameters {
  matrix[K, D] beta;
  matrix[K, J_age] alpha_age;
  matrix[K, J_eth] alpha_eth;
  matrix[K, J_edu] alpha_edu;
  vector[K] alpha;
  for (k in 1:(K - 1)) {
    alpha[k] = alpha_raw[k];
    alpha_age[k] = sigma_inter_eqn[1] * sigma_age[k] * eta_age[k];
    alpha_eth[k] = sigma_inter_eqn[2] * sigma_eth[k] * eta_eth[k];
    alpha_edu[k] = sigma_inter_eqn[3] * sigma_edu[k] * eta_edu[k];
  }
  alpha[K] = 0;
  beta = append_row(beta_raw, zeros);
  alpha_age[K] = zeros_age;
  alpha_eth[K] = zeros_eth;
  alpha_edu[K] = zeros_edu;
}
model {
  to_vector(beta_raw) ~ normal(0, 1);
  to_vector(eta_age) ~ normal(0, 1);
  to_vector(eta_eth) ~ normal(0, 1);
  to_vector(eta_edu) ~ normal(0, 1);
  sigma_age ~ normal(0, 1);
  sigma_eth ~ normal(0, 1);
  sigma_edu ~ normal(0, 1);
  alpha_raw ~ normal(0, 1);
  sigma_inter_eqn ~ normal(0, 1);
  {
    matrix[K, N] mu_logit;
    mu_logit = beta * t_X;
    for (n in 1:N)
      y[n] ~ categorical_logit(col(mu_logit,n)
                               + alpha + col(alpha_age,idx_age[n])
                               + col(alpha_eth,idx_eth[n]) + col(alpha_edu,idx_edu[n]));
  }
}