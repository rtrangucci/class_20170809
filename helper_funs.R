## Generate the cholesky factor of an LKJ
## correlation matrix
chol_lkj <- function(d, nu){
  shape1 <- rep(NA_real_, d - 1)
  shape2 <- rep(NA_real_, d - 1)
  alpha <- nu + (d - 2) / 2
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

brk_finder <- function(x) {
  x_range <- range(x)
  n <- nclass.Sturges(x)
  brks <- pretty(x_range, n = n, min.n = 1)
  return(brks)
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

intervals_of_interest <- function(ps, group_facs, cell_counts, ps_reps, probs) { 
  
  cell_N <- ps[,cell_counts] %>% data.matrix()
  ps$cell_N <- cell_N
  gpd_ps <- ps %>% ungroup() %>% group_by_(.dots = group_facs)
  gp_ind <- gpd_ps %>% group_indices()
  gpd_ps$gp_ind <- gp_ind
  gp_sums <- summarise_(gpd_ps, tot_pop = sum(cell_N), gp_ind = first(gp_ind))
  gp_nms <- select_(gp_sums, .dots = group_facs)
  gp_nms <- apply(gp_nms, 1, function(x) paste(x, collapse = '_'))
  
  gp_sums$lo = NA_real_
  gp_sums$med = NA_real_
  gp_sums$hi = NA_real_
  gp_num <- sort(unique(gp_ind))
  gp_map <- data.frame(nm = gp_nms, 
                       num = gp_num, stringsAsFactors = F)
  n_gps <- nrow(gp_map)
  ps_reps <- sweep(x = ps_reps, MARGIN = 2, STATS = ps$cell_N, FUN = '*')
  
  intervals <- as.data.frame(matrix(NA, nrow = n_gps, ncol = length(probs) + 1))
  dists <- list()
  for (gp_i in seq_along(gp_num)) {
    gp <- gp_map[gp_i,]
    gp_n <- gp$num
    gp_nm <- gp$nm
    sel <- which(gpd_ps$gp_ind == gp_n)
    weight_tot <- sum(ps$cell_N[sel])
    sub_vote <- ps_reps[, sel]
    vote_vec <- rowSums(sub_vote)/weight_tot
    intervals[gp_n,] <- round(c(quantile(vote_vec, probs = probs),mean(vote_vec)),3)
    dists[[gp_nm]] <- vote_vec
  }
  colnames(intervals) <- c(paste0(round(100*probs),'%'), "mean")
  intervals$group <- gp_map$nm
  return(list(intervals = intervals, 
              dists = dists)) 
}

softmax <- function(x) {exp(x) / sum(exp(x))}