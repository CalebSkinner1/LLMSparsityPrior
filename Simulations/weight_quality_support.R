# R script supplements weight_quality_sims.R to aid sampling

# load libraries ------------------------------------------------------------

suppressPackageStartupMessages({
  library("tidyverse")
  library("hypergeo")
  library("glmnet")
  library("parallel")
  library("furrr")
})

# generate weights -------------------------------------------------------

l1_weight_agreement <- function(true_gamma, weights) {
  # scale weights
  scaled_weights <- (weights - min(weights)) / (max(weights) - min(weights))

  # l1 distance
  l1 <- mean(abs(true_gamma - scaled_weights))

  1 - l1
}

l2_weight_agreement <- function(true_gamma, weights) {
  # scale weights
  scaled_weights <- (weights - min(weights)) / (max(weights) - min(weights))

  # l2 distance
  l2 <- mean((true_gamma - scaled_weights)^2)

  1 - l2
}

pairwise_weight_agreement <- function(true_gamma, weights) {
  # ensure true_gamma and weights have the same length
  if (length(true_gamma) != length(weights)) {
    stop("true_gamma and weights must have the same number of elements")
  }

  true_gamma_diff_mat <- outer(true_gamma, true_gamma, FUN = "-")
  true_gamma_mat <- ifelse(true_gamma_diff_mat > 0, 1, 0)

  weights_diff_mat <- outer(weights, weights, FUN = "-")
  weights_mat <- ifelse(weights_diff_mat > 0, 1, 0)

  # compute difference of true_gamma pairwise matrix and weights_diff pairwise matrix
  diff_mat <- abs(true_gamma_mat - weights_mat) # 1 if agreement differs, 0 if same
  total_disagreement <- sum(colSums(diff_mat))
  # compute the total number of pairwise interactions, ignoring the off diagonal and dismissing half the matrix
  total_cells <- length(true_gamma) * (length(true_gamma) - 1) / 2
  percent_agreement <- (total_cells - total_disagreement) / total_cells
  percent_agreement
}

ROC_weight_agreement <- function(true_gamma, weights){
  scaled_weights <- (weights - min(weights)) / (max(weights) - min(weights))

  ROC <- pROC::roc(true_gamma, scaled_weights, levels = c("0", "1"), direction = "<")
  auc_value <- pROC::auc(ROC)

  auc_value[1]
}

# generate weights (takes in l1 distance weight agreement and returns weights)
generate_weights <- function(
  phi, # weight agreement for l1 distance
  true_gamma
) {
  if (phi == 1) {
    weight_prop <- c(1, 0, 0, 0)
  } else {
    mu <- 3 * (1 - phi)
    # solve for ratio (geometric difference between weight proportions)
    f <- function(r) {
      (3 - mu) * r^3 + (2 - mu) * r^2 + (1 - mu) * r - mu
    }

    solution <- uniroot(f, interval = c(0, 1))
    r <- solution$root
    c <- mu / (r + 2 * r^2 + 3 * r^3)

    # proportion of weights that are (correct, 1 off, 2 off, 3 off)
    weight_prop <- c(c, c * r, c * r^2, c * r^3)
  }

  # noise (inactive) covariates
  p_0 <- sum(true_gamma == 0)
  weight_dist_0 <- weight_prop * p_0
  counts_0 <- floor(weight_dist_0)
  remainder_0 <- p_0 - sum(counts_0)
  if (remainder_0 > 0) {
    decimals_0 <- weight_dist_0 - counts_0
    # find indices of the highest decimals
    extra_indices_0 <- order(decimals_0, decreasing = TRUE)[1:remainder_0]
    counts_0[extra_indices_0] <- counts_0[extra_indices_0] + 1
  }

  weights_0 <- rep(
    1 + seq_along(counts_0),
    times = round(counts_0)
  )

  # signal (active) covariates
  p_1 <- sum(true_gamma == 1)
  weight_dist_1 <- weight_prop * p_1
  counts_1 <- floor(weight_dist_1)
  remainder_1 <- p_1 - sum(counts_1)
  if (remainder_1 > 0) {
    decimals_1 <- weight_dist_1 - counts_1
    # find indices of the highest decimals
    extra_indices_1 <- order(decimals_1, decreasing = TRUE)[1:remainder_1]
    counts_1[extra_indices_1] <- counts_1[extra_indices_1] + 1
  }

  weights_1 <- rep(
    6 - seq_along(counts_1),
    times = round(counts_1)
  )
  weights <- numeric(length(true_gamma))
  weights[which(true_gamma == 0)] <- weights_0
  weights[which(true_gamma == 1)] <- weights_1

  weights
}

# helper functions for simulation --------------------------------------------------------

compile_model_metrics <- function(model_object, weights, alpha, beta) {
  signal <- beta != 0

  if (class(model_object) == "numeric") {
    # i.e. lasso objects
    gamma_predict <- as.vector(model_object != 0)[-1]
    l1 <- sum(abs(as.vector(model_object) - c(alpha, beta))) # l1 loss of beta
  } else {
    # mcmc objects
    gamma_predict <- model_object$gamma
    l1 <- sum(abs(model_object$beta - c(alpha, beta))) # l1 loss of beta
  }

  # compute metrics
  fp <- length(which(gamma_predict[!signal] > 0.5)) # false positives of MPM
  fn <- length(which(gamma_predict[signal] < 0.5)) # false negatives
  tp <- sum(signal) - fn # true positives
  tn <- sum(!signal) - fp # true negatives
  precision <- if_else(tp == 0, 0, tp / (tp + fp)) # precision
  recall <- if_else(tp == 0, 0, tp / (tp + fn)) # recall

  # store metrics
  tibble(gamma_predict, weights, signal) %>%
    group_by(weights, signal) %>%
    summarize(coef_mean = mean(gamma_predict), .groups = "keep") %>%
    arrange(signal) %>%
    mutate(group_type = str_c("w", weights, "s", as.numeric(signal))) %>%
    ungroup() %>%
    select(group_type, coef_mean) %>%
    pivot_wider(names_from = group_type, values_from = coef_mean) %>%
    mutate(
      "f1" = if_else(
        precision + recall == 0,
        0,
        2 * (precision * recall) / (precision + recall)
      ),
      "l1" = l1,
      "fp" = fp,
      "fn" = fn
    )
}

# Load MCMC Samplers ----------------------------------------------------------------
# source("MCMC Samplers/LSP_regression_fixed_s.R")
# source("MCMC Samplers/LSP_regression_random_s.R")
lsp_fixed_log_posterior <- function(
  Z,
  Z_gram,
  y,
  tau,
  gamma,
  a_sigma,
  b_sigma,
  theta,
  n
) {
  n_gam <- ncol(Z) # cardinality of selected design matrix

  Q <- Z_gram + diag(n_gam) / tau
  cholQ <- chol(Q)
  log_detQ <- 2 * sum(log(diag(cholQ)))

  # compute model prior
  model_prior <- sum(gamma * log(theta) + (1 - gamma) * log(1 - theta))

  y_Z <- crossprod(Z, y)
  quadratic_term <- t(y_Z) %*% chol2inv(cholQ) %*% y_Z

  n_gam / 2 * log(tau) -
    .5 * log_detQ -
    (n / 2 + a_sigma) * log(sum(y^2) - quadratic_term + b_sigma) +
    model_prior
}

# compute log acceptance rate: log(p(gamma_new|data)) - log(p(gamma_old|data))
lsp_fixed_log_acceptance_rate <- function(
  Z_old,
  Z_old_gram,
  Z_new,
  Z_new_gram,
  y,
  gamma_new,
  gamma_old,
  tau,
  a_sigma,
  b_sigma,
  theta,
  n
) {
  lsp_fixed_log_posterior(
    Z_new,
    Z_new_gram,
    y,
    tau,
    gamma_new,
    a_sigma,
    b_sigma,
    theta,
    n
  ) -
    lsp_fixed_log_posterior(
      Z_old,
      Z_old_gram,
      y,
      tau,
      gamma_old,
      a_sigma,
      b_sigma,
      theta,
      n
    )
}

# run discrete spike and slab sampler for regression. With confidence = 0, reduces to traditional spike and slab.
# With confidence = 1, prior sparsity is entirely determined by the weights
lsp_fixed_gibbs_sampler <- function(
  X,
  y,
  weights = NULL,
  c = NA,
  eta = 0,
  sparsity,
  a_sigma,
  b_sigma,
  tau,
  iter = 10000,
  burn_in = 5000,
  thin = 1,
  prob_add = 1 / 3,
  prob_delete = 1 / 3,
  init_weights = TRUE,
  return_samples = TRUE
) {
  if (is.null(weights)) {
    c <- 0 # if no weights, then assign zero confidence and zero eta
    eta <- 0
  }

  p <- ncol(X)
  n <- nrow(X)

  # number of values of eta and c in grid
  K <- length(c)*length(eta)

  # all c, eta values
  cross_eta_c <- expand.grid(eta = eta, c = c)

  # create space for theta_mat
  theta_mat <- matrix(0, nrow = nrow(cross_eta_c), ncol = p)

  # eta and c are fixed
  if (length(c) == 1 & length(eta) == 1) {
    if(eta == 0 | c == 0){ 
      init_weights <- FALSE
      theta_mat[1,] <- rep(sparsity, p)
    }else{
      # create vector for prior model probability
      raw_theta <- sparsity * c * (weights^eta) / mean(weights^eta) + (1 - c) * sparsity

      theta_mat[1,] <- pmin(pmax(raw_theta, 1e-4), 1 - 1e-4)
    }
  } else { # eta and c are random
      for(k in 1:nrow(cross_eta_c)){
        c_k <- cross_eta_c$c[k]
        eta_k <- cross_eta_c$eta[k]

        raw_theta <- sparsity * c_k * (weights^eta_k) /mean((weights^eta_k)) + (1 - c_k) * sparsity

        # constrain theta to be less than 1
        capped_theta <- pmin(pmax(raw_theta, 1e-4), 1 - 1e-4)

        theta_mat[k, ] <- capped_theta        
      }
    }

  # number of iter left after thinning/burn_in
  n_keep <- ceiling((iter - burn_in) / thin)

  if(return_samples){
    # create space for gamma, beta, invsigma^2, acc
    gam_store <- matrix(0, nrow = n_keep, ncol = p)
    beta_store <- matrix(0, nrow = n_keep, ncol = p + 1)
    invsigma_2_store <- rep(0, n_keep)
    acc_store <- rep(0, n_keep)
    eta_store <- rep(0, n_keep)
    c_store <- rep(0, n_keep)
  } else{
    # reduction for memory: only store means
    gam_mean <- rep(0, p)
    beta_mean <- rep(0, p + 1)
    invsigma_2_mean <- 0
    acc_mean <- 0
    eta_mean <- 0
    c_mean <- 0
  }

  # generate initial values of gamma (in a smart way)
  gam_current <- rep(0, p)

  if (init_weights) {
    # initialize by taking largest p*sparsity weights from LLM
    gam_current[order(weights, decreasing = TRUE)[
      1:max(2, ceiling(sparsity * p))
    ]] <- 1
  } else {
    # initialize by taking largest p*sparsity correlations with y (marginal correlation screening)
    gam_current[order(abs(cor(X, y)), decreasing = TRUE)[
      1:max(2, ceiling(sparsity * p))
    ]] <- 1
  }

  fit <- glmnet::glmnet(
    X[, which(gam_current == 1)],
    y,
    alpha = 0,
    lambda = 1
  )
  beta_current <- as.vector(coef(fit))

  invsigma_2_current <- 1 /
    (1 /
      n *
      sum((y - cbind(1, X[, which(gam_current == 1)]) %*% beta_current)^2))
  
  # find current index of eta and c in discrete uniform grid
  if(K > 1){
    c_eta_idx_current <- sample(K, 1)
  }else{
    c_eta_idx_current <- 1
  }
  
  # begin iterations
  for (i in 1:iter) {
    # propose a candidate gamma
    gam_prop <- gam_current
    selected_gam <- which(gam_prop == 1)
    removed_gam <- which(gam_prop == 0)

    Z_old <- cbind(1, X[, selected_gam])
    Z_old_gram <- crossprod(Z_old)

    # employ ADS random search to edit gam_prop
    unif_gam <- runif(1)
    if (length(selected_gam) == 0) {
      unif_gam <- 0.5
    } # ensure that if none are selected, we will add

    log_prop_ratio <- 0
    current_model_size <- sum(gam_prop)

    # with prob_delete, randomly remove one gamma
    if (unif_gam < prob_delete || length(removed_gam) == 0) {
      chosen <- sample(selected_gam)[1]
      gam_prop[chosen] <- 0
      log_prop_ratio <- log(prob_add) - log(prob_delete) + log(current_model_size) -
        log(p - current_model_size + 1)
    } else if (unif_gam < prob_delete + prob_add) {
      # with prob_add, randomly add one gamma
      chosen <- sample(removed_gam)[1]
      gam_prop[chosen] <- 1
      log_prop_ratio <- log(prob_delete) - log(prob_add) + log(p - current_model_size) -
        log(current_model_size + 1)
    } else {
      # else swap
      chosen1 <- sample(removed_gam)[1]
      chosen2 <- sample(selected_gam)[1]
      gam_prop[chosen1] <- 1
      gam_prop[chosen2] <- 0
    }
    selected_gam <- which(gam_prop == 1)

    Z_new <- cbind(1, X[, selected_gam])
    Z_new_gram <- crossprod(Z_new)

    # compute log acceptance rate
    logacc <- lsp_fixed_log_acceptance_rate(
      Z_old,
      Z_old_gram,
      Z_new,
      Z_new_gram,
      y,
      gam_prop,
      gam_current,
      tau,
      a_sigma,
      b_sigma,
      theta_mat[c_eta_idx_current,],
      n
    ) +
      log_prop_ratio
    if (log(runif(1)) < logacc[[1]]) {
      # insert gamma draw
      gam_current <- gam_prop

      # draw beta vector
      chol_mat <- chol(
        invsigma_2_current *
          Z_new_gram +
          (invsigma_2_current / tau) * diag(ncol(Z_new_gram))
      )
      invQ <- chol2inv(chol_mat)
      l <- invsigma_2_current * crossprod(Z_new, y)
      beta_gamma <- MASS::mvrnorm(1, invQ %*% l, invQ)
      beta_current <- rep(0, p + 1)
      beta_current[c(1, selected_gam + 1)] <- beta_gamma

      # draw invsigma_2
      invsigma_2_current <- rgamma(
        1,
        shape = n / 2 + a_sigma,
        rate = 1 / 2 * sum((y - Z_new %*% beta_gamma)^2) + b_sigma
      )

      # count acceptances
      acc <- 1
    } else {
      # gam_current is maintained

      # draw beta vector
      chol_mat <- chol(
        invsigma_2_current *
          Z_old_gram +
          (invsigma_2_current / tau) * diag(ncol(Z_old_gram))
      )
      invQ <- chol2inv(chol_mat)
      l <- invsigma_2_current * crossprod(Z_old, y)
      beta_gamma <- MASS::mvrnorm(1, invQ %*% l, invQ)
      beta_current <- rep(0, p + 1)
      beta_current[c(1, which(gam_current == 1) + 1)] <- beta_gamma

      # draw invsigma_2
      invsigma_2_current <- rgamma(
        1,
        shape = n / 2 + a_sigma,
        rate = 1 / 2 * sum((y - Z_old %*% beta_gamma)^2) + b_sigma
      )

      # count acceptances
      acc <- 0
    }

    # draw new c and eta based on gamma

    # unnormalized log probability for all K states of eta and c
    W <- numeric(K)
    for(k in 1:K){
      W[k] <- sum(gam_current*log(theta_mat[k,]) + (1-gam_current)*log(1 - theta_mat[k,]))
    }
    # normalize probabilities
    # log-sum-exp trick to prevent NaN underflow
    pi_eta_c <- exp(W - max(W))/sum(exp(W - max(W)))
    c_eta_idx_current <- sample(K, 1, prob = pi_eta_c)

    # store parameters
    if (i > burn_in && (i - burn_in) %% thin == 0) {
      if(return_samples){
        store_i <- (i - burn_in) / thin # index

        gam_store[store_i, ] <- gam_current
        beta_store[store_i, ] <- beta_current
        invsigma_2_store[store_i] <- invsigma_2_current
        eta_store[store_i] <- cross_eta_c$eta[c_eta_idx_current]
        c_store[store_i] <- cross_eta_c$c[c_eta_idx_current]
        acc_store[store_i] <- acc
      }else{
        gam_mean <- gam_mean + gam_current/n_keep
        beta_mean <- beta_mean + beta_current/n_keep
        invsigma_2_mean <- invsigma_2_mean + invsigma_2_current/n_keep
        eta_mean <- eta_mean + cross_eta_c$eta[c_eta_idx_current]/n_keep
        c_mean <- c_mean + cross_eta_c$c[c_eta_idx_current]/n_keep
        acc_mean <- acc_mean + acc/n_keep
      }
    }
  }

  if(return_samples){
    list(
      "beta" = beta_store,
      "gamma" = gam_store,
      "invsigma_2" = invsigma_2_store,
      "eta" = eta_store,
      "c" = c_store,
      "accs" = acc_store)
  }else{
    list(
      "beta" = beta_mean,
      "gamma" = gam_mean,
      "invsigma_2" = invsigma_2_mean,
      "eta" = eta_mean,
      "c" = c_mean,
      "accs" = acc_mean)
  }
}

compute_log_prior_gamma <- function(gamma, s, u) {
  theta <- s * u # compute theta vector

  # compute log probability of gamma conditional on s
  log_prob <- sum(gamma * log(theta) + (1 - gamma) * log(1 - theta))

  log_prob
}

# compute unnormalized log-posterior density (marginalizing out beta and sigma)
# (Z is (1, X)), does not include some constants that will be cancelled in log_acceptance_rate
lsp_random_log_posterior <- function(
  Z,
  Z_gram,
  y,
  tau,
  gamma,
  a_sigma,
  b_sigma,
  s,
  u,
  n
) {
  n_gam <- ncol(Z) # cardinality of selected design matrix

  Q <- Z_gram + diag(n_gam) / tau
  cholQ <- chol(Q)
  log_detQ <- 2 * sum(log(diag(cholQ)))

  model_prior <- compute_log_prior_gamma(gamma, s, u)

  y_Z <- crossprod(Z, y)
  quadratic_term <- t(y_Z) %*% chol2inv(cholQ) %*% y_Z

  n_gam /2 *log(tau) - .5 * log_detQ -
    (n / 2 + a_sigma) *
      log(sum(y^2) - quadratic_term + b_sigma) +
    model_prior
}

# compute log acceptance rate: log(p(gamma_new|data)) - log(p(gamma_old|data))
lsp_random_log_acceptance_rate <- function(
  Z_old,
  Z_old_gram,
  Z_new,
  Z_new_gram,
  y,
  gamma_new,
  gamma_old,
  tau,
  a_sigma,
  b_sigma,
  s,
  u,
  n
) {
  lsp_random_log_posterior(
    Z_new,
    Z_new_gram,
    y,
    tau,
    gamma_new,
    a_sigma,
    b_sigma,
    s,
    u,
    n
  ) -
    lsp_random_log_posterior(
      Z_old,
      Z_old_gram,
      y,
      tau,
      gamma_old,
      a_sigma,
      b_sigma,
      s,
      u,
      n
    )
}

lsp_random_gibbs_sampler <- function(
  X,
  y,
  weights = NULL,
  c = NA,
  eta = 0,
  a_sigma,
  b_sigma,
  tau,
  a_s = 1,
  b_s = NA,
  s_proposal_sigma = 1,
  iter = 10000,
  burn_in = 5000,
  thin = 1,
  prob_add = 1 / 3,
  prob_delete = 1 / 3,
  init_weights = TRUE,
  return_samples = TRUE
) {
  if (is.null(weights)) {
    c <- 0 # if no weights, then clearly there is no confidence in them
    eta <- 0
  }

  p <- ncol(X)
  n <- nrow(X)

  # number of values of eta and c in grid
  K <- length(c)*length(eta)

  # all c, eta values
  cross_eta_c <- expand.grid(eta = eta, c = c)

  # create space for u_mat (u*sparsity = theta)
  u_mat <- matrix(0, nrow = nrow(cross_eta_c), ncol = p)

  # eta and c are fixed
  if (length(c) == 1 & length(eta) == 1) {
    if(eta == 0 | c == 0){ 
      init_weights <- FALSE
      u_mat[1,] <- rep(1, p)
    }else{
      # create vector for prior model probability
      u_mat[1,] <- c * (weights^eta) / mean(weights^eta) + (1 - c)
    }
  } else { # eta and c are random
      for(k in 1:nrow(cross_eta_c)){
        c_k <- cross_eta_c$c[k]
        eta_k <- cross_eta_c$eta[k]

        u_mat[k, ] <- c_k * (weights^eta_k) /mean((weights^eta_k)) + (1 - c_k)
      }
    }

  # per recommendation of rockova-george
  if (is.na(b_s)) {
    b_s <- p
  }

  # number of models left after thinning/burn_in
  n_keep <- ceiling((iter - burn_in) / thin)

  if(return_samples){
    # create space for gamma, beta, invsigma^2, acc
    gam_store <- matrix(0, nrow = n_keep, ncol = p)
    beta_store <- matrix(0, nrow = n_keep, ncol = p + 1)
    invsigma_2_store <- rep(0, n_keep)
    acc_store <- rep(0, n_keep)
    eta_store <- rep(0, n_keep)
    c_store <- rep(0, n_keep)
    s_store <- rep(0, n_keep)
    acc_s_store <- rep(0, n_keep)
  } else{
    # reduction for memory: only store means
    gam_mean <- rep(0, p)
    beta_mean <- rep(0, p + 1)
    invsigma_2_mean <- 0
    acc_mean <- 0
    eta_mean <- 0
    c_mean <- 0
    s_mean <- 0
    acc_s_mean <- 0
  }

  # generate initial values of gamma, s, sigma^2 (in a smart way)
  gam_current <- rep(0, p)

  s_current <- a_s / (a_s + b_s)

  if (init_weights) {
    # initialize by taking largest p*prior_mean_s weights from LLM
    gam_current[order(weights, decreasing = TRUE)[
      1:max(2, ceiling(s_current * p))
    ]] <- 1
  } else {
    # initialize by taking largest p*prior_mean_s correlations with y (marginal correlation screening)
    gam_current[order(abs(cor(X, y)), decreasing = TRUE)[
      1:max(2, ceiling(s_current * p))
    ]] <- 1
  }

  fit <- glmnet::glmnet(
    X[, which(gam_current == 1)],
    y,
    alpha = 0,
    lambda = 1
  )
  beta_current <- as.vector(coef(fit))

  invsigma_2_current <- 1 /
    (1 /
      n *
      sum((y - cbind(1, X[, which(gam_current == 1)]) %*% beta_current)^2))
  
  # find current index of eta and c in discrete uniform grid
  if(K > 1){
    c_eta_idx_current <- sample(K, 1)
  }else{
    c_eta_idx_current <- 1
  }
  
  # begin iterations
  for (i in 1:iter) {
    # propose a candidate gamma
    gam_prop <- gam_current
    selected_gam <- which(gam_prop == 1)
    removed_gam <- which(gam_prop == 0)

    Z_old <- cbind(1, X[, selected_gam])
    Z_old_gram <- crossprod(Z_old)

    # employ ADS random search to edit gam_prop
    unif_gam <- runif(1)
    if (length(selected_gam) == 0) {
      unif_gam <- 0.5
    } # ensure that if none are selected, we will add

    log_prop_ratio <- 0
    current_model_size <- sum(gam_prop)

    # with prob_delete, randomly remove one gamma
    if (unif_gam < prob_delete || length(removed_gam) == 0) {
      chosen <- sample(selected_gam)[1]
      gam_prop[chosen] <- 0
      log_prop_ratio <- log(prob_add) - log(prob_delete) + log(current_model_size) -
        log(p - current_model_size + 1)
    } else if (unif_gam < prob_delete + prob_add) {
      # with prob_add, randomly add one gamma
      chosen <- sample(removed_gam)[1]
      gam_prop[chosen] <- 1
      log_prop_ratio <-  log(prob_delete) - log(prob_add) +
        log(p - current_model_size) - log(current_model_size + 1)
    } else {
      # else swap
      chosen1 <- sample(removed_gam)[1]
      chosen2 <- sample(selected_gam)[1]
      gam_prop[chosen1] <- 1
      gam_prop[chosen2] <- 0
    }
    selected_gam <- which(gam_prop == 1)

    Z_new <- cbind(1, X[, selected_gam])
    Z_new_gram <- crossprod(Z_new)

    # compute log acceptance rate
    logacc <- lsp_random_log_acceptance_rate(
      Z_old,
      Z_old_gram,
      Z_new,
      Z_new_gram,
      y,
      gam_prop,
      gam_current,
      tau,
      a_sigma,
      b_sigma,
      s_current,
      u_mat[c_eta_idx_current,],
      n
    ) +
      log_prop_ratio
    if (log(runif(1)) < logacc[[1]]) {
      # insert gamma draw
      gam_current <- gam_prop

      # draw beta vector
      chol_mat <- chol(
        invsigma_2_current *
          Z_new_gram +
          (invsigma_2_current / tau) * diag(ncol(Z_new_gram))
      )
      invQ <- chol2inv(chol_mat)
      l <- invsigma_2_current * crossprod(Z_new, y)
      beta_gamma <- MASS::mvrnorm(1, invQ %*% l, invQ)
      beta_current <- rep(0, p + 1)
      beta_current[c(1, selected_gam + 1)] <- beta_gamma

      # draw invsigma_2
      invsigma_2_current <- rgamma(
        1,
        shape = n / 2 + a_sigma,
        rate = 1 / 2 * sum((y - Z_new %*% beta_gamma)^2) + b_sigma
      )

      # count acceptances
      acc <- 1
    } else {
      # gam_current is maintained

      # draw beta vector
      chol_mat <- chol(
        invsigma_2_current *
          Z_old_gram +
          (invsigma_2_current / tau) * diag(ncol(Z_old_gram))
      )
      invQ <- chol2inv(chol_mat)
      l <- invsigma_2_current * crossprod(Z_old, y)
      beta_gamma <- MASS::mvrnorm(1, invQ %*% l, invQ)
      beta_current <- rep(0, p + 1)
      beta_current[c(1, which(gam_current == 1) + 1)] <- beta_gamma

      # draw invsigma_2
      invsigma_2_current <- rgamma(
        1,
        shape = n / 2 + a_sigma,
        rate = 1 / 2 * sum((y - Z_old %*% beta_gamma)^2) + b_sigma
      )

      # count acceptances
      acc <- 0
    }

    # draw new s via random walk on logit scale
    logit_s <- log(s_current / (1 - s_current))
    logit_s_new <- rnorm(1, mean = logit_s, sd = s_proposal_sigma)
    s_new <- 1 / (1 + exp(-logit_s_new))

    if (max(s_new * u_mat[c_eta_idx_current,]) >= 1) {
      accept_s <- FALSE
    } else {
      # compute posterior ratio for s
      log_prior_s_new <- dbeta(s_new, a_s, b_s, log = TRUE)
      log_prior_s_old <- dbeta(s_current, a_s, b_s, log = TRUE)

      log_lik_s_new <- compute_log_prior_gamma(gam_current, s_new, u_mat[c_eta_idx_current,])
      log_lik_s_old <- compute_log_prior_gamma(gam_current, s_current, u_mat[c_eta_idx_current,])

      # jacobian adjustment
      log_jacobian_new <- log(s_new) + log(1 - s_new)
      log_jacobian_old <- log(s_current) + log(1 - s_current)

      # total Metrpolis-Hastings ratio for s proposal
      log_ratio_s <- (log_lik_s_new + log_prior_s_new + log_jacobian_new) -
        (log_lik_s_old + log_prior_s_old + log_jacobian_old)

      if (log(runif(1)) < log_ratio_s) {
        s_current <- s_new
        accept_s <- TRUE
      } else {
        accept_s <- FALSE
      }
    }

    # draw new c and eta based on gamma

    # unnormalized log probability for all K states of eta and c
    W <- numeric(K)
    for(k in 1:K){
      theta_k <- pmin(pmax(s_current*u_mat[k,], 1e-4), 1-1e-4)
      W[k] <- sum(gam_current*log(theta_k) + (1-gam_current)*log(1 - theta_k))
    }
    # normalize probabilities
    # log-sum-exp trick to prevent NaN underflow
    pi_eta_c <- exp(W - max(W))/sum(exp(W - max(W)))
    c_eta_idx_current <- sample(K, 1, prob = pi_eta_c)

    # store parameters
    if (i > burn_in && (i - burn_in) %% thin == 0) {
      if(return_samples){
        store_i <- (i - burn_in) / thin # index

        gam_store[store_i, ] <- gam_current
        beta_store[store_i, ] <- beta_current
        invsigma_2_store[store_i] <- invsigma_2_current
        eta_store[store_i] <- cross_eta_c$eta[c_eta_idx_current]
        c_store[store_i] <- cross_eta_c$c[c_eta_idx_current]
        s_store[store_i] <- s_current
        acc_store[store_i] <- acc
        acc_s_store[store_i] <- accept_s
      }else{
        gam_mean <- gam_mean + gam_current/n_keep
        beta_mean <- beta_mean + beta_current/n_keep
        invsigma_2_mean <- invsigma_2_mean + invsigma_2_current/n_keep
        eta_mean <- eta_mean + cross_eta_c$eta[c_eta_idx_current]/n_keep
        c_mean <- c_mean + cross_eta_c$c[c_eta_idx_current]/n_keep
        s_mean <- s_mean + s_current/n_keep
        acc_mean <- acc_mean + acc/n_keep
        acc_s_mean <- acc_s_mean + accept_s/n_keep
      }
    }
  }

  if(return_samples){
    list(
      "beta" = beta_store,
      "gamma" = gam_store,
      "invsigma_2" = invsigma_2_store,
      "eta" = eta_store,
      "c" = c_store,
      "s" = s_store,
      "accs" = acc_store,
      "acc_s" = acc_s_store)
  }else{
    list(
      "beta" = beta_mean,
      "gamma" = gam_mean,
      "invsigma_2" = invsigma_2_mean,
      "eta" = eta_mean,
      "c" = c_mean,
      "s" = s_mean,
      "accs" = acc_mean,
      "acc_s" = acc_s_mean)
  }
}

# LLM-Lasso ---------------------------------------------------------------
# code may be found at https://github.com/pilancilab/LLM-Lasso, lightly edited for functionality

# Scale X_test with X_train's center/scale; guard zero-variance cols
.scale_like_train <- function(X_train, X_new = NULL) {
  X_train <- as.matrix(X_train)
  center <- colMeans(X_train)
  scalev <- apply(X_train, 2, sd)
  scalev[!is.finite(scalev) | scalev == 0] <- 1
  X_train_sc <- scale(X_train, center = center, scale = scalev)
  X_new_sc <- if (!is.null(X_new)) {
    scale(as.matrix(X_new), center = center, scale = scalev)
  } else {
    NULL
  }
  list(X_train = X_train_sc, X_new = X_new_sc, center = center, scale = scalev)
}

# Align (optionally named) weights to X's column order; validate positivity
.align_and_check_weights <- function(weights, X) {
  if (!is.null(names(weights))) {
    if (!all(colnames(X) %in% names(weights))) {
      missing <- setdiff(colnames(X), names(weights))
      stop(
        "weights are named, but missing entries for features: ",
        paste(missing, collapse = ", ")
      )
    }
    w <- as.numeric(weights[colnames(X)])
  } else {
    w <- as.numeric(weights)
    if (length(w) != ncol(X)) stop("length(weights) must equal ncol(X)")
  }
  if (any(!is.finite(w)) || any(w <= 0)) {
    stop("All weights must be positive and finite")
  }
  pmax(w, 1e-8)
}

# Area between CV curves across matched sparsity (larger = better)
cve <- function(cvm, non_zero, ref_cvm, ref_non_zero) {
  df1 <- tibble(ref_non_zero, ref_cvm) %>%
    group_by(ref_non_zero) %>%
    summarise(ref_cvm = min(ref_cvm), .groups = 'drop') %>%
    arrange(ref_non_zero)
  df2 <- tibble(non_zero, cvm) %>%
    group_by(non_zero) %>%
    summarise(cvm = min(cvm), .groups = 'drop') %>%
    arrange(non_zero)

  interp <- stats::approx(
    x = df1$ref_non_zero,
    y = df1$ref_cvm,
    xout = df2$non_zero,
    method = "linear",
    rule = 2
  )

  n <- length(df2$non_zero)
  if (n < 2) {
    return(0)
  }

  area <- 0
  for (i in 1:(n - 1)) {
    width <- df2$non_zero[[i + 1]] - df2$non_zero[[i]]
    height <- ((interp$y[[i]] - df2$cvm[[i]]) +
      (interp$y[[i + 1]] - df2$cvm[[i + 1]])) /
      2
    area <- area + width * height
  }
  area
}

llm_lasso_simp <- function(
  X_train,
  y_train,
  weights,
  folds_cv = 5,
  elastic_net = 1,
  max_imp_pow = 10,
  lambda_min_ratio = 0.01,
  regression = TRUE,
  multinomial = FALSE,
  type_measure = NULL,
  use_lambda_1se = FALSE
) {
  glm_family <-
    if (multinomial) {
      "multinomial"
    } else if (regression) {
      "gaussian"
    } else {
      "binomial"
    }
  if (is.null(type_measure)) {
    type_measure <- if (glm_family == "gaussian") "mse" else "class"
  }
  if (
    glm_family %in%
      c("binomial", "multinomial") &&
      !type_measure %in% c("class", "deviance")
  ) {
    stop('For classification, type_measure must be "class" or "deviance".')
  }

  # Align labels for classification
  if (glm_family %in% c("binomial", "multinomial")) {
    y_train <- if (is.factor(y_train)) y_train else factor(y_train)
  }

  X_train_sc <- .scale_like_train(X_train)$X_train
  w <- .align_and_check_weights(weights, X_train_sc)

  pf_list <- lapply(0:max_imp_pow, function(i) 1 / (w^i))
  pf_names <- paste0("1/imp^", 0:max_imp_pow)

  ref_cvm <- ref_nz <- NULL
  best_area <- -Inf
  best_name <- NULL
  best_pf <- NULL

  for (k in seq_along(pf_list)) {
    cv <- glmnet::cv.glmnet(
      x = X_train_sc,
      y = y_train,
      family = glm_family,
      alpha = elastic_net,
      penalty.factor = pf_list[[k]],
      nfolds = folds_cv,
      lambda.min.ratio = lambda_min_ratio,
      standardize = FALSE,
      type.measure = type_measure
    )
    if (is.null(ref_cvm)) {
      ref_cvm <- cv$cvm
    }
    if (is.null(ref_nz)) {
      ref_nz <- cv$nzero
    }
    a <- cve(cv$cvm, cv$nzero, ref_cvm, ref_nz)
    if (a > best_area) {
      best_area <- a
      best_name <- pf_names[k]
      best_pf <- pf_list[[k]]
    }
  }

  cv_best <- glmnet::cv.glmnet(
    x = X_train_sc,
    y = y_train,
    family = glm_family,
    alpha = elastic_net,
    penalty.factor = best_pf,
    nfolds = folds_cv,
    lambda.min.ratio = lambda_min_ratio,
    standardize = FALSE,
    type.measure = type_measure
  )
  s_choice <- if (use_lambda_1se) "lambda.1se" else "lambda.min"
  lam <- if (use_lambda_1se) cv_best$lambda.1se else cv_best$lambda.min

  if (glm_family != "multinomial") {
    co <- as.numeric(coef(cv_best, s = s_choice))
    n_features <- sum(co[-1] != 0)
    coef_obj <- co
  } else {
    co_list <- coef(cv_best, s = s_choice) # list (one coef matrix per class)
    nz_by_class <- lapply(co_list, function(cm) as.numeric(cm[-1, 1] != 0))
    feat_nonzero <- Reduce("|", nz_by_class)
    n_features <- sum(feat_nonzero)
    coef_obj <- co_list
  }

  list(
    model = best_name,
    algo = "LLM-Lasso",
    method = lam,
    coef = coef_obj, # numeric vector (gaussian/binomial) or list per class (multinomial)
    n_features = n_features
  )
}

# simulation function -----------------------------------------------------
# this function runs the simulation for baseline methods
baseline_data_sim_function <- function(seed, n){
  set.seed(seed)
  # generate X and coefficients
  X <- MASS::mvrnorm(n, mu = rep(0, p), Xvar * cov_mat)
  beta <- c(rep(0, p - s), rep(effect_size, s))
  alpha <- effect_size

  # generate y
  y <- X %*% beta + alpha + rnorm(n, 0, sd = y_sd)

  # run lasso
  lasso_results <- as.vector(coef(glmnet::glmnet(
    X,
    y,
    alpha = 1,
    lambda = glmnet::cv.glmnet(X, y, alpha = 1)$lambda.min
  )))

    # horseshoe prior
  hs_fit <- Mhorseshoe::approx_horseshoe(
    y = y, 
    X = cbind(1, X), 
    burn = 10000, 
    iter = 5000
)
  hs_coef <- hs_fit$BetaHat
  rm(hs_fit)
  gc()
  
  # Spike-and-Slab Lasso
  sslasso_fit <- SSLASSO::SSLASSO(X = X, y = y)
  model_index <- floor(ncol(sslasso_fit$beta)/2)
  sslasso_intercept <- sslasso_fit$intercept[model_index]
  sslasso_beta <- sslasso_fit$beta[, model_index]
  sslasso_coefs <- c(sslasso_intercept, sslasso_beta)

  rm(sslasso_fit)
  gc()

  baseline_fits <- list(
    "lasso" = lasso_results,
    "horseshoe" = hs_coef,
    "ss_lasso" = sslasso_coefs)
  
  if (fixed_s == TRUE) {
    
    # run standard discrete spike and slab
    baseline_fits$`ss, fixed s` <- lsp_fixed_gibbs_sampler(
      X,
      y,
      c = 0,
      eta = 0,
      sparsity = sparsity,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau,
      iter = iter,
      burn_in = burn_in,
      init_weights = FALSE,
      return_samples = FALSE
    )
  }
  # can also evaluate LSP with random sparsity
  if (random_s == TRUE) {
    # run standard discrete spike and slab with random s
    baseline_fits$`ss, random s` <- lsp_random_gibbs_sampler(
      X,
      y,
      c = 0,
      eta = 0,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau,
      iter = iter,
      burn_in = burn_in,
      init_weights = FALSE,
      return_samples = FALSE
    )
  }

  # return generated data and baseline models
  list(
    data = list(X = X, y = y, alpha = alpha, beta = beta),
    baselines = baseline_fits
  )
}


sim_function <- function(baseline_fits, weights) {
  # generate X and coefficients
  X <- baseline_fits$data$X
  y <- baseline_fits$data$y
  beta <- baseline_fits$data$beta
  alpha <- baseline_fits$data$alpha

  all_fits <- baseline_fits$baselines

  # run llm-lasso
  all_fits$`llm-lasso` <- llm_lasso_simp(
    X_train = X,
    y_train = y,
    weights = weights,
    elastic_net = 1,
    regression = TRUE
  )$coef

  if (fixed_s == TRUE) {
    
    # run LSP with fixed sparsity
    all_fits$`lsp, fixed s` <- lsp_fixed_gibbs_sampler(
      X,
      y,
      weights,
      sparsity = sparsity,
      c = confidence_range,
      eta = eta_range,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau,
      iter = iter,
      burn_in = burn_in,
      init_weights = TRUE,
      return_samples = FALSE
    )
  }

  # can also evaluate LSP with random sparsity
  if (random_s == TRUE) {

    # run LSP with random sparsity
    all_fits$`lsp, random s` <- lsp_random_gibbs_sampler(
      X,
      y,
      weights,
      c = confidence_range,
      eta = eta_range,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau,
      iter = iter,
      burn_in = burn_in,
      init_weights = TRUE,
      return_samples = FALSE
    )
  }

  # return metrics
  metrics <- map(
    all_fits,
    ~ compile_model_metrics(.x, weights, alpha, beta)
  )

  return(metrics)
}
