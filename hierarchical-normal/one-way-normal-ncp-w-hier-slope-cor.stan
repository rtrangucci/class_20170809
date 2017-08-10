data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1, upper=J> idx_J[N];
  vector[N] y;
  vector[N] X;
}
parameters {
  matrix[2,J] eta;
  real<lower=0> sigma_y;
  vector<lower=0>[2] sigma;
  real alpha;
  real gamma;
  cholesky_factor_corr[2] L_Omega;
}
transformed parameters {
  vector[J] theta;
  vector[J] beta;
  
  {
    matrix[J,2] dummy = (diag_pre_multiply(sigma, L_Omega) * eta)';
    theta = col(dummy,1);
    beta = gamma + col(dummy,2);
  }
}
model {
  y ~ normal(alpha + theta[idx_J] + X .* beta[idx_J], sigma_y);
  
  sigma_y ~ cauchy(0, 10);
  alpha ~ normal(100, 20);
  to_vector(eta) ~ normal(0, 1); 
  gamma ~ normal(0, 10);
  sigma ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(3.0);
}
generated quantities {
  matrix[2, 2] Omega = L_Omega * L_Omega';
}