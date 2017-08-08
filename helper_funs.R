## Generate the cholesky factor of an LKJ
## correlation matrix
chol_lkj <- function(d, eta){
  shape1 <- rep(NA_real_, d - 1)
  shape2 <- rep(NA_real_, d - 1)
  alpha <- eta + (d - 2) / 2
  shape1[1] <- alpha
  shape2[1] <- alpha
  for(i in 2:(d - 1)){
    alpha <- alpha - 1 / 2
    shape2[i] <- alpha
    shape1[i] <- i / 2
  }
  r2 <- rbeta(d - 1, shape1, shape2)
  L <- matrix(0,d,d)
  L[1,1] <- 1
  L[2,1] <- 2 * r2[1] - 1
  L[2,2] <- sqrt(1 - L[2,1] ^ 2)
  for(m in 2:(d - 1)){
    l <- rnorm(m)
    scale <- sqrt(r2[m] / (t(l) %*% l)[1,])
    L[m + 1,1:m] <- l * scale
    L[m + 1, m + 1] <- sqrt(1 - r2[m])
  }
  return(L)
}

n_divergent <- function(fit) {
  sampler_diagnostics <- get_sampler_params(fit, inc_warmup = FALSE)
  return(sapply(sampler_diagnostics, function(x) sum(x[,'divergent__'])))
}

fac2int <- function(cat_var) {
  return(as.numeric(as.factor(cat_var)))
}

n_levels <- function(cat_var) {
  return(length(unique(cat_var)))
}