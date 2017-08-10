data { 
  int<lower=2> K;
  int<lower=0> N;
  int<lower=0> y[N, K];
  int<lower=1> J_age;
  int<lower=1> J_state;
  int<lower=1> G;
  int<lower=1, upper=J_age> idx_age[N];
  int<lower=1, upper=J_state> idx_state[N];
}
transformed data {
  row_vector[J_age] zeros_age = rep_row_vector(0, J_age);
  row_vector[J_state] zeros_state = rep_row_vector(0, J_state);
}
parameters {
  vector[K - 1] alpha_raw;
  matrix[K - 1, J_age] eta_age;
  matrix[K - 1, J_state] eta_state;
  vector<lower=0>[K - 1] sigma_age;
  vector<lower=0>[K - 1] sigma_state;
  vector<lower=0>[G] sigma_inter_eqn;
}
transformed parameters {
  matrix[K, J_age] alpha_age;
  matrix[K, J_state] alpha_state;
  vector[K] alpha;
  for (k in 1:(K - 1)) {
    alpha[k] = alpha_raw[k];
    alpha_age[k] = sigma_inter_eqn[1] * sigma_age[k] * eta_age[k];
    alpha_state[k] = sigma_inter_eqn[2] * sigma_state[k] * eta_state[k];
  }
  alpha[K] = 0;
  alpha_age[K] = zeros_age;
  alpha_state[K] = zeros_state;
}
model {
  to_vector(eta_age) ~ normal(0, 1);
  to_vector(eta_state) ~ normal(0, 1);
  sigma_inter_eqn ~ normal(0, 1);
  sigma_age ~ normal(0, 1);
  sigma_state ~ normal(0, 1);
  alpha_raw ~ normal(0, 1);
  for (n in 1:N)
    y[n] ~ multinomial(softmax(col(alpha_age,idx_age[n])
                               + col(alpha_state,idx_state[n])));
}
