data { 
  int<lower=2> K;
  int<lower=0> N;
  int<lower=1> D;
  int<lower=1,upper=K> y[N];
  matrix[N, D] X;
}
transformed data {
  row_vector[D] zeros = rep_row_vector(0, D);
  matrix[D, N] t_X = X';
}
parameters {
  matrix[K - 1, D] beta_raw;
  vector[K - 1] alpha_raw;
}
transformed parameters {
  matrix[K, D] beta;
  vector[K] alpha;
  beta = append_row(beta_raw, zeros);
  for (k in 1:(K - 1))
    alpha[k] = alpha_raw[k];
  alpha[K] = 0;
}
model {
  to_vector(beta_raw) ~ normal(0, 1);
  alpha_raw ~ normal(0, 1);
  {
    matrix[K, N] mu_logit;
    mu_logit = beta * t_X;
    for (n in 1:N)
      y[n] ~ categorical_logit(col(mu_logit,n) + alpha);
  }
}