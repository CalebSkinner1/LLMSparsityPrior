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

ROC_weight_agreement <- function(true_gamma, weights) {
  scaled_weights <- (weights - min(weights)) / (max(weights) - min(weights))

  ROC <- pROC::roc(
    true_gamma,
    scaled_weights,
    levels = c("0", "1"),
    direction = "<"
  )
  auc_value <- pROC::auc(ROC)

  auc_value[1]
}

# generate weights (takes in l1 distance weight agreement and returns weights)
generate_weights <- function(
  phi, # weight agreement for l1 distance
  true_gamma,
  categories = 5
) {
  if (phi == 1) {
    weight_prop <- c(1, rep(0, categories - 1))
  } else {
    mu <- (categories - 1) * (1 - phi)
    # solve for ratio (geometric difference between weight proportions)
    f <- function(r) {
      k <- 1:(categories - 1)
      sum((k - mu) * r^k) - mu
    }

    solution <- uniroot(f, interval = c(0, 1))
    r <- solution$root

    k <- 1:(categories - 1)
    c <- mu / sum(k * r^k)

    # proportion of weights that are (correct, 1 off, 2 off, 3 off)
    weight_prop <- c * (r^(0:(categories - 1)))
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
    seq_along(counts_0),
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
    (categories + 1) - seq_along(counts_1),
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
lsp_fixed_ss_log_posterior <- function(
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
  quadratic_term <- as.numeric(t(y_Z) %*% chol2inv(cholQ) %*% y_Z)

  n_gam /
    2 *
    log(tau) -
    .5 * log_detQ -
    (n / 2 + a_sigma) * log(.5 * (sum(y^2) - quadratic_term) + b_sigma) +
    model_prior
}

# compute log acceptance rate: log(p(gamma_new|data)) - log(p(gamma_old|data))
lsp_fixed_ss_log_acceptance_rate <- function(
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
  lsp_fixed_ss_log_posterior(
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
    lsp_fixed_ss_log_posterior(
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
lsp_fixed_ss_gibbs_sampler <- function(
  X,
  y,
  weights = NULL,
  c = NA,
  eta = NULL,
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
  } else if (is.null(eta)) {
    eta <- 0
    step_size <- 1
    # Calculate initial bound to ensure it starts < 1
    theta_bound <- s * max(weights)^eta / mean(weights^eta)

    while (eta <= 20) {
      eta <- eta + step_size # step forward
      theta_bound <- s * max(weights)^eta / mean(weights^eta)

      # if threshold is crossed, backtrack with smaller steps
      if (theta_bound >= 1) {
        # Step back to the last safe value
        eta <- eta - step_size

        # Decrease the step size for finer searching
        if (step_size == 1) {
          step_size <- 0.1
        } else if (step_size == 0.1) {
          step_size <- 0.01
        } else {
          break
        }
      }
    }
    # generate eta space
    eta <- seq(1, eta, length.out = 10)
  }

  p <- ncol(X)
  n <- nrow(X)

  # number of values of eta and c in grid
  K <- length(c) * length(eta)

  # all c, eta values
  cross_eta_c <- expand.grid(eta = eta, c = c)

  # create space for theta_mat
  theta_mat <- matrix(0, nrow = nrow(cross_eta_c), ncol = p)

  # eta and c are fixed
  if (length(c) == 1 & length(eta) == 1) {
    if (eta == 0 | c == 0) {
      init_weights <- FALSE
      if (length(sparsity) == 1) {
        theta_mat[1, ] <- rep(sparsity, p)
      } else {
        # insert inclusion probabilities directly
        theta_mat[1, ] <- sparsity
      }
    } else {
      # create vector for prior model probability
      raw_theta <- sparsity *
        c *
        (weights^eta) /
        mean(weights^eta) +
        (1 - c) * sparsity

      theta_mat[1, ] <- pmin(pmax(raw_theta, 1e-4), 1 - 1e-4)
    }
  } else {
    # eta and c are random
    for (k in 1:nrow(cross_eta_c)) {
      c_k <- cross_eta_c$c[k]
      eta_k <- cross_eta_c$eta[k]

      raw_theta <- sparsity *
        c_k *
        (weights^eta_k) /
        mean((weights^eta_k)) +
        (1 - c_k) * sparsity

      # constrain theta to be less than 1
      capped_theta <- pmin(pmax(raw_theta, 1e-4), 1 - 1e-4)

      theta_mat[k, ] <- capped_theta
    }
  }

  # number of iter left after thinning/burn_in
  n_keep <- ceiling((iter - burn_in) / thin)

  if (return_samples) {
    # create space for gamma, beta, invsigma^2, acc
    gam_store <- matrix(0, nrow = n_keep, ncol = p)
    beta_store <- matrix(0, nrow = n_keep, ncol = p + 1)
    invsigma_2_store <- rep(0, n_keep)
    acc_store <- rep(0, n_keep)
    eta_store <- rep(0, n_keep)
    c_store <- rep(0, n_keep)
  } else {
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
  if (K > 1) {
    c_eta_idx_current <- sample(K, 1)
  } else {
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
      log_prop_ratio <- log(prob_add) -
        log(prob_delete) +
        log(current_model_size) -
        log(p - current_model_size + 1)
    } else if (unif_gam < prob_delete + prob_add) {
      # with prob_add, randomly add one gamma
      chosen <- sample(removed_gam)[1]
      gam_prop[chosen] <- 1
      log_prop_ratio <- log(prob_delete) -
        log(prob_add) +
        log(p - current_model_size) -
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
    logacc <- lsp_fixed_ss_log_acceptance_rate(
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
      theta_mat[c_eta_idx_current, ],
      n
    ) +
      log_prop_ratio
    if (log(runif(1)) < logacc) {
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
    for (k in 1:K) {
      W[k] <- sum(
        gam_current *
          log(theta_mat[k, ]) +
          (1 - gam_current) * log(1 - theta_mat[k, ])
      )
    }
    # normalize probabilities
    # log-sum-exp trick to prevent NaN underflow
    pi_eta_c <- exp(W - max(W)) / sum(exp(W - max(W)))
    c_eta_idx_current <- sample(K, 1, prob = pi_eta_c)

    # store parameters
    if (i > burn_in && (i - burn_in) %% thin == 0) {
      c_val <- cross_eta_c$c[c_eta_idx_current]
      if (return_samples) {
        store_i <- (i - burn_in) / thin # index

        gam_store[store_i, ] <- gam_current
        beta_store[store_i, ] <- beta_current
        invsigma_2_store[store_i] <- invsigma_2_current
        c_store[store_i] <- c_val
        eta_store[store_i] <- if (c_val == 0) {
          0
        } else {
          cross_eta_c$eta[c_eta_idx_current]
        }
        acc_store[store_i] <- acc
      } else {
        gam_mean <- gam_mean + gam_current / n_keep
        beta_mean <- beta_mean + beta_current / n_keep
        invsigma_2_mean <- invsigma_2_mean + invsigma_2_current / n_keep
        c_mean <- c_mean + c_val / n_keep
        eta_mean <- if (c_val == 0) {
          eta_mean
        } else {
          eta_mean + cross_eta_c$eta[c_eta_idx_current] / n_keep
        }
        acc_mean <- acc_mean + acc / n_keep
      }
    }
  }

  if (return_samples) {
    list(
      "beta" = beta_store,
      "gamma" = gam_store,
      "invsigma_2" = invsigma_2_store,
      "eta" = eta_store,
      "c" = c_store,
      "accs" = acc_store
    )
  } else {
    list(
      "beta" = beta_mean,
      "gamma" = gam_mean,
      "invsigma_2" = invsigma_2_mean,
      "eta" = eta_mean,
      "c" = c_mean,
      "accs" = acc_mean
    )
  }
}

# compute model prior conditional on s
compute_log_prior_gamma <- function(gamma, s, u) {
  theta <- pmax(1e-10, pmin(s * u, 1 - 1e-10)) # compute theta vector

  # compute log probability of gamma conditional on s
  log_prob <- sum(gamma * log(theta) + (1 - gamma) * log(1 - theta))

  log_prob
}

# compute unnormalized log-posterior density (marginalizing out beta and sigma)
# (Z is (1, X)), does not include some constants that will be cancelled in log_acceptance_rate
lsp_random_ss_log_posterior <- function(
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
  quadratic_term <- as.numeric(t(y_Z) %*% chol2inv(cholQ) %*% y_Z)

  n_gam /
    2 *
    log(tau) -
    .5 * log_detQ -
    (n / 2 + a_sigma) * log(.5 * (sum(y^2) - quadratic_term) + b_sigma) +
    model_prior
}

# compute log acceptance rate: log(p(gamma_new|data)) - log(p(gamma_old|data))
lsp_random_ss_log_acceptance_rate <- function(
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
  lsp_random_ss_log_posterior(
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
    lsp_random_ss_log_posterior(
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

lsp_random_ss_gibbs_sampler <- function(
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
    c <- 0 # if no weights, then assign zero confidence and zero eta
    eta <- 0
  } else if (is.null(eta)) {
    # if missing eta space, assign by default
    eta <- 0
    step_size <- 1
    # Calculate initial bound to ensure it starts < 1; hard code sparsity to be 0.05 for upper bound
    theta_bound <- 0.05 * max(weights)^eta / mean(weights^eta)

    while (eta <= 19) {
      eta <- eta + step_size # step forward
      theta_bound <- 0.05 * max(weights)^eta / mean(weights^eta)

      # if threshold is crossed, backtrack with smaller steps
      if (theta_bound >= 1) {
        # Step back to the last safe value
        eta <- eta - step_size

        # Decrease the step size for finer searching
        if (step_size == 1) {
          step_size <- 0.1
        } else if (step_size == 0.1) {
          step_size <- 0.01
        } else {
          break
        }
      }
    }
    # generate eta space
    eta <- seq(1, eta, length.out = 10)
  }

  p <- ncol(X)
  n <- nrow(X)

  # number of values of eta and c in grid
  K <- length(c) * length(eta)

  # all c, eta values
  cross_eta_c <- expand.grid(eta = eta, c = c)

  # create space for u_mat (u*sparsity = theta)
  u_mat <- matrix(0, nrow = nrow(cross_eta_c), ncol = p)

  # eta and c are fixed
  if (length(c) == 1 & length(eta) == 1) {
    if (eta == 0 | c == 0) {
      init_weights <- FALSE
      u_mat[1, ] <- rep(1, p)
    } else {
      # create vector for prior model probability
      u_mat[1, ] <- c * (weights^eta) / mean(weights^eta) + (1 - c)
    }
  } else {
    # eta and c are random
    for (k in 1:nrow(cross_eta_c)) {
      c_k <- cross_eta_c$c[k]
      eta_k <- cross_eta_c$eta[k]

      u_mat[k, ] <- c_k * (weights^eta_k) / mean((weights^eta_k)) + (1 - c_k)
    }
  }

  # per recommendation of rockova-george
  if (is.na(b_s)) {
    b_s <- p
  }

  # number of models left after thinning/burn_in
  n_keep <- ceiling((iter - burn_in) / thin)

  if (return_samples) {
    # create space for gamma, beta, invsigma^2, acc
    gam_store <- matrix(0, nrow = n_keep, ncol = p)
    beta_store <- matrix(0, nrow = n_keep, ncol = p + 1)
    invsigma_2_store <- rep(0, n_keep)
    acc_store <- rep(0, n_keep)
    eta_store <- rep(0, n_keep)
    c_store <- rep(0, n_keep)
    s_store <- rep(0, n_keep)
    acc_s_store <- rep(0, n_keep)
  } else {
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
  if (K > 1) {
    c_eta_idx_current <- sample(K, 1)
  } else {
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
      log_prop_ratio <- log(prob_add) -
        log(prob_delete) +
        log(current_model_size) -
        log(p - current_model_size + 1)
    } else if (unif_gam < prob_delete + prob_add) {
      # with prob_add, randomly add one gamma
      chosen <- sample(removed_gam)[1]
      gam_prop[chosen] <- 1
      log_prop_ratio <- log(prob_delete) -
        log(prob_add) +
        log(p - current_model_size) -
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
    logacc <- lsp_random_ss_log_acceptance_rate(
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
      u_mat[c_eta_idx_current, ],
      n
    ) +
      log_prop_ratio
    if (log(runif(1)) < logacc) {
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

    if (max(s_new * u_mat[c_eta_idx_current, ]) >= 1) {
      accept_s <- FALSE
    } else {
      # compute posterior ratio for s
      log_prior_s_new <- dbeta(s_new, a_s, b_s, log = TRUE)
      log_prior_s_old <- dbeta(s_current, a_s, b_s, log = TRUE)

      log_lik_s_new <- compute_log_prior_gamma(
        gam_current,
        s_new,
        u_mat[c_eta_idx_current, ]
      )
      log_lik_s_old <- compute_log_prior_gamma(
        gam_current,
        s_current,
        u_mat[c_eta_idx_current, ]
      )

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
    for (k in 1:K) {
      theta_k <- pmin(pmax(s_current * u_mat[k, ], 1e-4), 1 - 1e-4)
      W[k] <- sum(
        gam_current * log(theta_k) + (1 - gam_current) * log(1 - theta_k)
      )
    }
    # normalize probabilities
    # log-sum-exp trick to prevent NaN underflow
    pi_eta_c <- exp(W - max(W)) / sum(exp(W - max(W)))
    c_eta_idx_current <- sample(K, 1, prob = pi_eta_c)

    # store parameters
    if (i > burn_in && (i - burn_in) %% thin == 0) {
      c_val <- cross_eta_c$c[c_eta_idx_current]
      if (return_samples) {
        store_i <- (i - burn_in) / thin # index

        gam_store[store_i, ] <- gam_current
        beta_store[store_i, ] <- beta_current
        invsigma_2_store[store_i] <- invsigma_2_current
        c_store[store_i] <- c_val
        eta_store[store_i] <- if (c_val == 0) {
          0
        } else {
          cross_eta_c$eta[c_eta_idx_current]
        }
        s_store[store_i] <- s_current
        acc_store[store_i] <- acc
        acc_s_store[store_i] <- accept_s
      } else {
        gam_mean <- gam_mean + gam_current / n_keep
        beta_mean <- beta_mean + beta_current / n_keep
        invsigma_2_mean <- invsigma_2_mean + invsigma_2_current / n_keep
        c_mean <- c_mean + c_val / n_keep
        eta_mean <- if (c_val == 0) {
          eta_mean
        } else {
          eta_mean + cross_eta_c$eta[c_eta_idx_current] / n_keep
        }
        s_mean <- s_mean + s_current / n_keep
        acc_mean <- acc_mean + acc / n_keep
        acc_s_mean <- acc_s_mean + accept_s / n_keep
      }
    }
  }

  if (return_samples) {
    list(
      "beta" = beta_store,
      "gamma" = gam_store,
      "invsigma_2" = invsigma_2_store,
      "eta" = eta_store,
      "c" = c_store,
      "s" = s_store,
      "accs" = acc_store,
      "acc_s" = acc_s_store
    )
  } else {
    list(
      "beta" = beta_mean,
      "gamma" = gam_mean,
      "invsigma_2" = invsigma_2_mean,
      "eta" = eta_mean,
      "c" = c_mean,
      "s" = s_mean,
      "accs" = acc_mean,
      "acc_s" = acc_s_mean
    )
  }
}

# Load Spike-and-Slab LASSO - fixed ----------------------------------------------
lsp_fixed_ssl_map <- function(
  X,
  y,
  weights,
  s_fixed = 0.5,
  E_space = c(1, 2, 3),
  c_space = c(0, 1),
  lambda1 = NULL,
  lambda0s = NULL,
  variance = "unknown",
  max_iter = 500,
  eps = 0.001
) {
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(weights)) {
    weights <- rep(1, p)
    E_space <- 0
    c_space <- 0
  }

  # ── Standardize X
  std_X <- standardize(X)
  XX <- std_X$XX
  mean_y <- mean(y)
  yy <- y - mean_y

  # default lambda1 is 1
  if (is.null(lambda1)) {
    lambda1 <- 1
  }
  # default lambda0 path is 1:n
  if (is.null(lambda0s)) {
    lambda0s <- seq(lambda1, n, length.out = 100)
  }

  # Sigma initialization
  df <- 3
  sigquant <- 0.9
  sigest <- sd(yy)
  qchi <- qchisq(1 - sigquant, df)
  ncp <- sigest^2 * qchi / df
  min_sigma2 <- sigest^2 / n
  sigma_init <- if (variance == "unknown") sqrt(df * ncp / (df + 2)) else 1.0

  initialbeta <- rep(0, p)
  count_max <- 10L

  # ── Fit ────────────────────────────────────────────────────────────────
  map_results <- find_MAP_fixed_hyperparams(
    X = XX,
    y = yy,
    weights = weights,
    s_fixed = s_fixed,
    E_space = E_space,
    c_space = c_space,
    lambda1 = lambda1,
    lambda0s = lambda0s,
    initialbeta = initialbeta,
    variance = variance,
    sigma = sigma_init,
    min_sigma2 = min_sigma2,
    eps = eps,
    max_iter = max_iter,
    count_max = count_max
  )

  # ── Post-processing ───────────
  for (l in seq_along(map_results)) {
    map_results[[l]]$intercept <- mean_y -
      as.numeric(crossprod(std_X$c, map_results[[l]]$beta))
    names(map_results[[l]]$beta) <- colnames(X)
  }

  return(map_results)
}

standardize <- function(X) {
  c_vec <- colMeans(X)
  XX <- sweep(X, 2, c_vec, "-") # centre only, matching scale(..., scale=FALSE)
  list(XX = XX, c = c_vec)
}

construct_v_vec <- function(weights, c_val, eta) {
  if (c_val == 0) {
    return(rep(1.0, length(weights)))
  }
  w_eta <- weights^eta
  c_val * (w_eta / mean(w_eta)) + (1 - c_val)
}

# Vectorised over a vector of theta_j values
threshold_func_vec <- function(theta_vec, sigma2, lambda1, lambda0, xnorm_vec) {
  p <- length(theta_vec)
  if (lambda1 == lambda0) {
    return(rep(sigma2 * lambda1, p))
  }

  # pstar at 0
  ratio0 <- ((1 - theta_vec) / theta_vec) * (lambda0 / lambda1)
  pstar0 <- 1 / (1 + ratio0)
  lstar0 <- pstar0 * lambda1 + (1 - pstar0) * lambda0

  # g_func at 0
  g0 <- (lstar0 - lambda1)^2 + 2 * xnorm_vec / sigma2 * log(pstar0)

  result <- ifelse(
    g0 > 0,
    sqrt(2 * xnorm_vec * sigma2 * log(1 / pstar0)) + sigma2 * lambda1,
    sigma2 * lstar0
  )
  result
}

update_sigma2 <- function(r) sum(r^2) / (length(r))

# Scalar SSL thresholding
SSL_thresholding <- function(
  z,
  beta_val,
  lambda0,
  lambda1,
  theta_j,
  norm,
  delta,
  sigma2
) {
  if (abs(z) <= delta) {
    return(0)
  }
  if (lambda1 == lambda0) {
    lambda <- lambda1
  } else {
    ratio <- ((1 - theta_j) / theta_j) *
      (lambda0 / lambda1) *
      exp(-abs(beta_val) * (lambda0 - lambda1))
    pstar_val <- 1 / (1 + ratio)
    lambda <- pstar_val * lambda1 + (1 - pstar_val) * lambda0
  }
  temp <- abs(z) - sigma2 * lambda
  if (temp > 0) sign(z) * temp / norm else 0
}

# log posterior for evaluating c, eta (fixed s)
compute_log_posterior_fixed <- function(
  y,
  X,
  beta_final,
  sigma_final,
  v_vec,
  lambda1,
  lambda0_final,
  s_fixed
) {
  n <- length(y)
  p <- length(beta_final)

  theta_v <- pmax(pmin(s_fixed * v_vec, 1 - 1e-4), 1e-4)

  residuals <- y - as.vector(X %*% beta_final)
  rss <- sum(residuals^2)
  sig2 <- sigma_final^2
  log_lik <- -(n / 2) * log(sig2) - rss / (2 * sig2)

  ab <- abs(beta_final)
  term1 <- theta_v * (lambda1 / 2) * exp(-lambda1 * ab)
  term2 <- (1 - theta_v) * (lambda0_final / 2) * exp(-lambda0_final * ab)
  log_prior_beta <- sum(log(pmax(term1 + term2, 1e-300)))

  # log_prior_s is omitted because s is fixed (it acts as a constant)

  log_lik + log_prior_beta
}

lsp_ssl_fixed_descent <- function(
  X,
  y,
  initialbeta,
  variance = "unknown",
  lambda1,
  lambda0s,
  v_vec,
  s_fixed,
  sigma,
  min_sigma2,
  eps,
  max_iter,
  count_max
) {
  n <- nrow(X)
  p <- ncol(X)
  L <- length(lambda0s)

  xnorm <- colSums(X^2)
  b_mat <- matrix(0, nrow = p, ncol = L)
  sigmas <- rep(NA_real_, L)
  loss <- rep(NA_real_, L)
  iter_vec <- rep(0L, L)

  a_vec <- initialbeta
  newa <- initialbeta
  a_old <- numeric(p)

  e1 <- as.integer(a_vec != 0)
  e2 <- as.integer(a_vec != 0)

  r <- y - as.vector(X %*% a_vec)
  z <- as.vector(crossprod(X, r))

  thres <- min(n, max_iter)
  if (p < thres) {
    XTY <- as.vector(crossprod(X, y))
    XTX <- crossprod(X)
  } else {
    XTY <- XTX <- NULL
  }

  delta <- numeric(p)
  sigma2 <- sigma^2
  sigma2_init <- sigma^2
  estimate_sigma <- FALSE

  # Theta is fixed throughout the path based on s_fixed
  theta_vec <- pmax(pmin(s_fixed * v_vec, 1 - 1e-10), 1e-10)

  for (l in 1:L) {
    lambda0 <- lambda0s[l]

    # ── Sigma update ────────────────────────────────────────────────────
    if (l > 1L && variance == "unknown") {
      if (iter_vec[l - 1L] < 100L) {
        estimate_sigma <- TRUE
        sigma2 <- update_sigma2(r)
        if (sigma2 < min_sigma2) {
          sigma2 <- sigma2_init
          estimate_sigma <- FALSE
        }
      } else {
        estimate_sigma <- FALSE
        if (iter_vec[l - 1L] == max_iter) sigma2 <- sigma2_init
      }
    }

    # Delta always updates on new lambda0
    delta <- threshold_func_vec(theta_vec, sigma2, lambda1, lambda0, xnorm)
    e2 <- pmax(e2, as.integer(abs(z) > delta))
    counter <- 0L

    # ══ Outer while ═══════════════════════════════════════════════════════
    while (iter_vec[l] < max_iter) {
      # ── Active-set optimisation ────────────────────────────────────────
      while (iter_vec[l] < max_iter) {
        iter_vec[l] <- iter_vec[l] + 1L
        a_old[] <- a_vec

        for (j in 1:p) {
          counter <- counter + 1L # always increment, for all j

          if (e1[j]) {
            z[j] <- if (p >= thres) {
              crossprod(X[, j], r)[[1L]] + xnorm[j] * a_vec[j]
            } else {
              XTY[j] - crossprod(XTX[, j], newa)[[1L]] + xnorm[j] * a_vec[j]
            }

            b_val <- SSL_thresholding(
              z[j],
              a_vec[j],
              lambda0,
              lambda1,
              theta_vec[j],
              xnorm[j],
              delta[j],
              sigma2
            )
            b_mat[j, l] <- b_val

            if (p >= thres) {
              shift <- b_val - a_vec[j]
              if (shift != 0) r <- r - shift * X[, j]
            } else {
              newa[j] <- b_val
            }

            if (b_val == 0) e1[j] <- 0L
          }

          # Periodic update: fires every count_max variable visits
          if (counter == count_max) {
            if (variance == "unknown" && estimate_sigma) {
              sigma2 <- update_sigma2(r)
              if (sigma2 < min_sigma2) {
                sigma2 <- sigma2_init
              }
              # Re-evaluate delta only if sigma changed
              delta <- threshold_func_vec(
                theta_vec,
                sigma2,
                lambda1,
                lambda0,
                xnorm
              )
            }
            counter <- 0L
          }
        }

        # Sync a_vec
        if (p >= thres) {
          a_vec <- b_mat[, l]
        } else {
          a_vec <- newa
          newa <- a_vec
        }

        # Convergence check
        check_idx <- which(e1 == 1L | a_old != 0)
        converged_active <- TRUE
        if (length(check_idx) > 0L) {
          denom <- abs(a_vec[check_idx])
          zero_new <- denom == 0
          denom[zero_new] <- abs(a_old[check_idx][zero_new])
          denom[denom == 0] <- .Machine$double.eps
          if (any(abs(a_vec[check_idx] - a_old[check_idx]) / denom > eps)) {
            converged_active <- FALSE
          }
        }
        if (converged_active) break
      } # end active-set while

      # ── Strong-set violation scan ──────────────────────────────────────
      violations_strong <- 0L
      counter <- 0L
      strong_cands <- which(e1 == 0L & e2 == 1L)

      for (j in strong_cands) {
        z[j] <- if (p >= thres) {
          crossprod(X[, j], r)[[1L]] + xnorm[j] * a_vec[j]
        } else {
          XTY[j] - crossprod(XTX[, j], a_vec)[[1L]] + xnorm[j] * a_vec[j]
        }

        b_val <- SSL_thresholding(
          z[j],
          a_vec[j],
          lambda0,
          lambda1,
          theta_vec[j],
          xnorm[j],
          delta[j],
          sigma2
        )
        b_mat[j, l] <- b_val

        if (b_val != 0) {
          e1[j] <- 1L
          e2[j] <- 1L
          if (p >= thres) {
            r <- r - b_val * X[, j]
          }
          a_vec[j] <- b_val
          violations_strong <- violations_strong + 1L
          counter <- counter + 1L
        }

        if (counter == count_max) {
          if (variance == "unknown" && estimate_sigma) {
            sigma2 <- update_sigma2(r)
            if (sigma2 < min_sigma2) {
              sigma2 <- sigma2_init
            }
            delta <- threshold_func_vec(
              theta_vec,
              sigma2,
              lambda1,
              lambda0,
              xnorm
            )
          }
          counter <- 0L
        }
      }

      if (violations_strong > 0L) {
        next
      }

      # ── Rest violation scan ────────────────────────────────────────────
      violations_rest <- 0L
      counter <- 0L
      rest_cands <- which(e2 == 0L)

      for (j in rest_cands) {
        z[j] <- if (p >= thres) {
          crossprod(X[, j], r)[[1L]] + xnorm[j] * a_vec[j]
        } else {
          XTY[j] - crossprod(XTX[, j], a_vec)[[1L]] + xnorm[j] * a_vec[j]
        }

        b_val <- SSL_thresholding(
          z[j],
          a_vec[j],
          lambda0,
          lambda1,
          theta_vec[j],
          xnorm[j],
          delta[j],
          sigma2
        )
        b_mat[j, l] <- b_val

        if (b_val != 0) {
          e1[j] <- 1L
          e2[j] <- 1L
          if (p >= thres) {
            r <- r - b_val * X[, j]
          }
          a_vec[j] <- b_val
          violations_rest <- violations_rest + 1L
          counter <- counter + 1L
        }

        if (counter == count_max) {
          if (variance == "unknown" && estimate_sigma) {
            sigma2 <- update_sigma2(r)
            if (sigma2 < min_sigma2) {
              sigma2 <- sigma2_init
            }
            delta <- threshold_func_vec(
              theta_vec,
              sigma2,
              lambda1,
              lambda0,
              xnorm
            )
          }
          counter <- 0L
        }
      }

      if (violations_rest > 0L) {
        next
      }

      # ── Finalise ──────────────────────────────────────────────────────
      if (p < thres) {
        r <- y - as.vector(X %*% a_vec)
      }
      b_mat[, l] <- a_vec
      loss[l] <- sum(r^2)
      sigmas[l] <- sqrt(sigma2)
      break
    } # end outer while

    # Fallback
    if (is.na(loss[l])) {
      if (p < thres) {
        r <- y - as.vector(X %*% a_vec)
      }
      b_mat[, l] <- a_vec
      loss[l] <- sum(r^2)
      sigmas[l] <- sqrt(sigma2)
    }
  } # end lambda path

  list(beta = b_mat, loss = loss, iter = iter_vec, sigmas = sigmas)
}

# ── Grid search wrapper ───────────────────────────────────────────────────────

find_MAP_fixed_hyperparams <- function(
  X,
  y,
  weights,
  s_fixed,
  E_space,
  c_space,
  lambda1,
  lambda0s,
  initialbeta,
  variance,
  sigma,
  min_sigma2,
  eps,
  max_iter,
  count_max
) {
  L <- length(lambda0s)

  # Pre-compute all (c, eta) paths
  all_paths <- vector("list", length(c_space) * length(E_space))
  combo_idx <- 1L
  for (c_val in c_space) {
    for (eta_val in E_space) {
      v_vec <- construct_v_vec(weights, c_val, eta_val)
      res <- lsp_ssl_fixed_descent(
        X,
        y,
        v_vec = v_vec,
        s_fixed = s_fixed,
        variance = variance,
        lambda1 = lambda1,
        lambda0s = lambda0s,
        initialbeta = initialbeta,
        sigma = sigma,
        min_sigma2 = min_sigma2,
        eps = eps,
        max_iter = max_iter,
        count_max = count_max
      )
      all_paths[[combo_idx]] <- list(
        c_val = c_val,
        eta_val = eta_val,
        v_vec = v_vec,
        res = res
      )
      combo_idx <- combo_idx + 1L
    }
  }

  # Cross-sectional MAP selection
  best_path_results <- vector("list", L)

  for (l in 1:L) {
    best_score <- -Inf
    best_params_l <- list(
      lambda0 = lambda0s[l],
      c = NA,
      eta = NA,
      s = s_fixed,
      beta = NULL,
      sigma = NA,
      score = NA
    )

    for (path in all_paths) {
      beta_l <- path$res$beta[, l]
      sigma_l <- path$res$sigmas[l]

      if (is.na(sigma_l) || sigma_l <= 0 || anyNA(beta_l)) {
        next
      }

      current_score <- compute_log_posterior_fixed(
        y = y,
        X = X,
        beta_final = beta_l,
        sigma_final = sigma_l,
        v_vec = path$v_vec,
        lambda1 = lambda1,
        lambda0_final = lambda0s[l],
        s_fixed = s_fixed
      )

      if (!is.finite(current_score)) {
        next
      }

      if (current_score > best_score) {
        best_score <- current_score
        best_params_l$c <- path$c_val
        best_params_l$eta <- path$eta_val
        best_params_l$beta <- beta_l
        best_params_l$sigma <- sigma_l
        best_params_l$score <- current_score
      }
    }

    if (is.infinite(best_score)) {
      warning(sprintf("No valid score for lambda0[%d] = %.4f", l, lambda0s[l]))
    }
    best_path_results[[l]] <- best_params_l
  }

  best_path_results
}

# Spike-and-Slab LASSO - random ------------------------------------------
lsp_random_ssl_map <- function(
  X,
  y,
  weights,
  a_s = 1,
  b_s = NA,
  E_space = c(1, 2, 3),
  c_space = c(0, 1),
  lambda1 = NULL,
  lambda0s = NULL,
  variance = "unknown",
  max_iter = 500,
  eps = 0.001
) {
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(weights)) {
    weights <- rep(1, p)
    E_space <- 0
    c_space <- 0
  }

  if (is.na(b_s)) {
    b_s <- p
  }

  # ── Standardize X
  std_X <- standardize(X)
  XX <- std_X$XX
  mean_y <- mean(y)
  yy <- y - mean_y

  # default lambda1 is 1
  if (is.null(lambda1)) {
    lambda1 <- 1
  }
  # default lambda0 path is 1:n
  if (is.null(lambda0s)) {
    lambda0s <- seq(lambda1, n, length.out = 100)
  }

  # Sigma initialization
  df <- 3
  sigquant <- 0.9
  sigest <- sd(yy)
  qchi <- qchisq(1 - sigquant, df)
  ncp <- sigest^2 * qchi / df
  min_sigma2 <- sigest^2 / n
  sigma_init <- if (variance == "unknown") sqrt(df * ncp / (df + 2)) else 1.0

  initialbeta <- rep(0, p)
  count_max <- 10L

  # ── Fit ────────────────────────────────────────────────────────────────
  map_results <- find_MAP_hyperparams(
    X = XX,
    y = yy,
    weights = weights,
    a_s = a_s,
    b_s = b_s,
    E_space = E_space,
    c_space = c_space,
    lambda1 = lambda1,
    lambda0s = lambda0s,
    initialbeta = initialbeta,
    variance = variance,
    sigma = sigma_init,
    min_sigma2 = min_sigma2,
    eps = eps,
    max_iter = max_iter,
    count_max = count_max
  )

  # ── Post-processing ───────────
  for (l in seq_along(map_results)) {
    map_results[[l]]$intercept <- mean_y -
      as.numeric(crossprod(std_X$c, map_results[[l]]$beta))
    names(map_results[[l]]$beta) <- colnames(X)
  }

  return(map_results)
}

# Helper functions -------------------------------------------------------

standardize <- function(X) {
  c_vec <- colMeans(X)
  XX <- sweep(X, 2, c_vec, "-") # centre only, matching scale(..., scale=FALSE)
  list(XX = XX, c = c_vec)
}

expectation_approx <- function(beta_vec, a_s, b_s) {
  p <- length(beta_vec)
  max(0, min(1, (sum(beta_vec != 0) + a_s) / (a_s + b_s + p)))
}

construct_v_vec <- function(weights, c_val, eta) {
  if (c_val == 0) {
    return(rep(1.0, length(weights)))
  }
  w_eta <- weights^eta
  c_val * (w_eta / mean(w_eta)) + (1 - c_val)
}

# Vectorised over a vector of theta_j values
threshold_func_vec <- function(theta_vec, sigma2, lambda1, lambda0, xnorm_vec) {
  p <- length(theta_vec)
  if (lambda1 == lambda0) {
    return(rep(sigma2 * lambda1, p))
  }

  # pstar at 0
  ratio0 <- ((1 - theta_vec) / theta_vec) * (lambda0 / lambda1)
  pstar0 <- 1 / (1 + ratio0)
  lstar0 <- pstar0 * lambda1 + (1 - pstar0) * lambda0

  # g_func at 0
  g0 <- (lstar0 - lambda1)^2 + 2 * xnorm_vec / sigma2 * log(pstar0)

  result <- ifelse(
    g0 > 0,
    sqrt(2 * xnorm_vec * sigma2 * log(1 / pstar0)) + sigma2 * lambda1,
    sigma2 * lstar0
  )
  result
}

update_sigma2 <- function(r) sum(r^2) / (length(r))

# Scalar SSL thresholding
SSL_thresholding <- function(
  z,
  beta_val,
  lambda0,
  lambda1,
  theta_j,
  norm,
  delta,
  sigma2
) {
  if (abs(z) <= delta) {
    return(0)
  }
  if (lambda1 == lambda0) {
    lambda <- lambda1
  } else {
    ratio <- ((1 - theta_j) / theta_j) *
      (lambda0 / lambda1) *
      exp(-abs(beta_val) * (lambda0 - lambda1))
    pstar_val <- 1 / (1 + ratio)
    lambda <- pstar_val * lambda1 + (1 - pstar_val) * lambda0
  }
  temp <- abs(z) - sigma2 * lambda
  if (temp > 0) sign(z) * temp / norm else 0
}

# log posterior for evaluating c, eta
compute_log_posterior <- function(
  y,
  X,
  beta_final,
  sigma_final,
  v_vec,
  lambda1,
  lambda0_final,
  a_s,
  b_s
) {
  n <- length(y)
  p <- length(beta_final)

  s_val <- max(0, min(1, (sum(beta_final != 0) + a_s) / (a_s + b_s + p)))
  theta_v <- pmax(pmin(s_val * v_vec, 1 - 1e-4), 1e-4)

  residuals <- y - as.vector(X %*% beta_final)
  rss <- sum(residuals^2)
  sig2 <- sigma_final^2
  log_lik <- -(n / 2) * log(sig2) - rss / (2 * sig2)

  ab <- abs(beta_final)
  term1 <- theta_v * (lambda1 / 2) * exp(-lambda1 * ab)
  term2 <- (1 - theta_v) * (lambda0_final / 2) * exp(-lambda0_final * ab)
  log_prior_beta <- sum(log(pmax(term1 + term2, 1e-300)))

  s_safe <- max(min(s_val, 1 - 1e-10), 1e-10)
  log_prior_s <- (a_s - 1) * log(s_safe) + (b_s - 1) * log(1 - s_safe)

  log_lik + log_prior_beta + log_prior_s
}

# Core coordinate descent ------------------------------------------------

lsp_ssl_random_descent <- function(
  X,
  y,
  initialbeta,
  variance = "unknown",
  lambda1,
  lambda0s,
  v_vec,
  a_s,
  b_s,
  sigma,
  min_sigma2,
  eps,
  max_iter,
  count_max
) {
  n <- nrow(X)
  p <- ncol(X)
  L <- length(lambda0s)

  xnorm <- colSums(X^2)
  b_mat <- matrix(0, nrow = p, ncol = L)
  sigmas <- rep(NA_real_, L)
  loss <- rep(NA_real_, L)
  iter_vec <- rep(0L, L)
  s_path <- rep(NA_real_, L)

  a_vec <- initialbeta
  newa <- initialbeta
  a_old <- numeric(p)

  e1 <- as.integer(a_vec != 0)
  e2 <- as.integer(a_vec != 0)

  r <- y - as.vector(X %*% a_vec)
  z <- as.vector(crossprod(X, r))

  thres <- min(n, max_iter)
  if (p < thres) {
    XTY <- as.vector(crossprod(X, y))
    XTX <- crossprod(X)
  } else {
    XTY <- XTX <- NULL
  }

  delta <- numeric(p)
  sigma2 <- sigma^2
  sigma2_init <- sigma^2
  estimate_sigma <- FALSE
  theta_vec <- pmax(pmin(0.5 * v_vec, 1 - 1e-10), 1e-10)

  for (l in 1:L) {
    lambda0 <- lambda0s[l]

    # ── Sigma update ────────────────────────────────────────────────────
    if (l > 1L && variance == "unknown") {
      if (iter_vec[l - 1L] < 100L) {
        estimate_sigma <- TRUE
        sigma2 <- update_sigma2(r)
        if (sigma2 < min_sigma2) {
          sigma2 <- sigma2_init
          estimate_sigma <- FALSE
        }
      } else {
        estimate_sigma <- FALSE
        if (iter_vec[l - 1L] == max_iter) sigma2 <- sigma2_init
      }
    }

    # ── Theta / delta warm-start ─────────────────────────────────────────
    if (l > 1L) {
      s_val <- expectation_approx(b_mat[, l - 1L], a_s, b_s)
      theta_vec <- pmax(pmin(s_val * v_vec, 1 - 1e-10), 1e-10)
    }

    delta <- threshold_func_vec(theta_vec, sigma2, lambda1, lambda0, xnorm)
    e2 <- pmax(e2, as.integer(abs(z) > delta))
    counter <- 0L

    # ══ Outer while ═══════════════════════════════════════════════════════
    while (iter_vec[l] < max_iter) {
      # ── Active-set optimisation ────────────────────────────────────────
      while (iter_vec[l] < max_iter) {
        iter_vec[l] <- iter_vec[l] + 1L
        a_old[] <- a_vec

        for (j in 1:p) {
          counter <- counter + 1L # always increment, for all j

          if (e1[j]) {
            z[j] <- if (p >= thres) {
              crossprod(X[, j], r)[[1L]] + xnorm[j] * a_vec[j]
            } else {
              XTY[j] - crossprod(XTX[, j], newa)[[1L]] + xnorm[j] * a_vec[j]
            }

            b_val <- SSL_thresholding(
              z[j],
              a_vec[j],
              lambda0,
              lambda1,
              theta_vec[j],
              xnorm[j],
              delta[j],
              sigma2
            )
            b_mat[j, l] <- b_val

            if (p >= thres) {
              shift <- b_val - a_vec[j]
              if (shift != 0) r <- r - shift * X[, j]
            } else {
              newa[j] <- b_val
            }

            if (b_val == 0) e1[j] <- 0L
          }

          # Periodic update: fires every count_max variable visits
          # regardless of active set membership — matches C++ cadence
          if (counter == count_max) {
            if (variance == "unknown" && estimate_sigma) {
              sigma2 <- update_sigma2(r)
              if (sigma2 < min_sigma2) sigma2 <- sigma2_init
            }
            s_val <- expectation_approx(b_mat[, l], a_s, b_s)
            theta_vec <- pmax(pmin(s_val * v_vec, 1 - 1e-10), 1e-10)
            delta <- threshold_func_vec(
              theta_vec,
              sigma2,
              lambda1,
              lambda0,
              xnorm
            )
            counter <- 0L
          }
        }

        # Sync a_vec
        if (p >= thres) {
          a_vec <- b_mat[, l]
        } else {
          a_vec <- newa
          newa <- a_vec
        }

        # Convergence check
        check_idx <- which(e1 == 1L | a_old != 0)
        converged_active <- TRUE
        if (length(check_idx) > 0L) {
          denom <- abs(a_vec[check_idx])
          zero_new <- denom == 0
          denom[zero_new] <- abs(a_old[check_idx][zero_new])
          denom[denom == 0] <- .Machine$double.eps
          if (any(abs(a_vec[check_idx] - a_old[check_idx]) / denom > eps)) {
            converged_active <- FALSE
          }
        }
        if (converged_active) break
      } # end active-set while

      # ── Strong-set violation scan ──────────────────────────────────────
      violations_strong <- 0L
      counter <- 0L
      strong_cands <- which(e1 == 0L & e2 == 1L)

      for (j in strong_cands) {
        z[j] <- if (p >= thres) {
          crossprod(X[, j], r)[[1L]] + xnorm[j] * a_vec[j]
        } else {
          XTY[j] - crossprod(XTX[, j], a_vec)[[1L]] + xnorm[j] * a_vec[j]
        }

        b_val <- SSL_thresholding(
          z[j],
          a_vec[j],
          lambda0,
          lambda1,
          theta_vec[j],
          xnorm[j],
          delta[j],
          sigma2
        )
        b_mat[j, l] <- b_val

        if (b_val != 0) {
          e1[j] <- 1L
          e2[j] <- 1L
          if (p >= thres) {
            r <- r - b_val * X[, j]
          }
          a_vec[j] <- b_val
          violations_strong <- violations_strong + 1L
          counter <- counter + 1L
        }

        if (counter == count_max) {
          if (variance == "unknown" && estimate_sigma) {
            sigma2 <- update_sigma2(r)
            if (sigma2 < min_sigma2) sigma2 <- sigma2_init
          }
          s_val <- expectation_approx(b_mat[, l], a_s, b_s)
          theta_vec <- pmax(pmin(s_val * v_vec, 1 - 1e-10), 1e-10)
          delta <- threshold_func_vec(
            theta_vec,
            sigma2,
            lambda1,
            lambda0,
            xnorm
          )
          counter <- 0L
        }
      }

      if (violations_strong > 0L) {
        next
      }

      # ── Rest violation scan ────────────────────────────────────────────
      violations_rest <- 0L
      counter <- 0L
      rest_cands <- which(e2 == 0L)

      for (j in rest_cands) {
        z[j] <- if (p >= thres) {
          crossprod(X[, j], r)[[1L]] + xnorm[j] * a_vec[j]
        } else {
          XTY[j] - crossprod(XTX[, j], a_vec)[[1L]] + xnorm[j] * a_vec[j]
        }

        b_val <- SSL_thresholding(
          z[j],
          a_vec[j],
          lambda0,
          lambda1,
          theta_vec[j],
          xnorm[j],
          delta[j],
          sigma2
        )
        b_mat[j, l] <- b_val

        if (b_val != 0) {
          e1[j] <- 1L
          e2[j] <- 1L
          if (p >= thres) {
            r <- r - b_val * X[, j]
          }
          a_vec[j] <- b_val
          violations_rest <- violations_rest + 1L
          counter <- counter + 1L
        }

        if (counter == count_max) {
          if (variance == "unknown" && estimate_sigma) {
            sigma2 <- update_sigma2(r)
            if (sigma2 < min_sigma2) sigma2 <- sigma2_init
          }
          s_val <- expectation_approx(b_mat[, l], a_s, b_s)
          theta_vec <- pmax(pmin(s_val * v_vec, 1 - 1e-10), 1e-10)
          delta <- threshold_func_vec(
            theta_vec,
            sigma2,
            lambda1,
            lambda0,
            xnorm
          )
          counter <- 0L
        }
      }

      if (violations_rest > 0L) {
        next
      }

      # ── Finalise ──────────────────────────────────────────────────────
      if (p < thres) {
        r <- y - as.vector(X %*% a_vec)
      }
      b_mat[, l] <- a_vec
      loss[l] <- sum(r^2)
      sigmas[l] <- sqrt(sigma2)
      s_path[l] <- expectation_approx(b_mat[, l], a_s, b_s)
      break
    } # end outer while

    # Fallback
    if (is.na(loss[l])) {
      if (p < thres) {
        r <- y - as.vector(X %*% a_vec)
      }
      b_mat[, l] <- a_vec
      loss[l] <- sum(r^2)
      sigmas[l] <- sqrt(sigma2)
      s_path[l] <- expectation_approx(b_mat[, l], a_s, b_s)
    }
  } # end lambda path

  list(
    beta = b_mat,
    loss = loss,
    iter = iter_vec,
    sigmas = sigmas,
    s_path = s_path
  )
}

# ── Grid search wrapper ───────────────────────────────────────────────────────

find_MAP_hyperparams <- function(
  X,
  y,
  weights,
  a_s,
  b_s,
  E_space,
  c_space,
  lambda1,
  lambda0s,
  initialbeta,
  variance,
  sigma,
  min_sigma2,
  eps,
  max_iter,
  count_max
) {
  L <- length(lambda0s)

  # Pre-compute all (c, eta) paths
  all_paths <- vector("list", length(c_space) * length(E_space))
  combo_idx <- 1L
  for (c_val in c_space) {
    for (eta_val in E_space) {
      v_vec <- construct_v_vec(weights, c_val, eta_val)
      res <- lsp_ssl_random_descent(
        X,
        y,
        v_vec = v_vec,
        a_s = a_s,
        b_s = b_s,
        variance = variance,
        lambda1 = lambda1,
        lambda0s = lambda0s,
        initialbeta = initialbeta,
        sigma = sigma,
        min_sigma2 = min_sigma2,
        eps = eps,
        max_iter = max_iter,
        count_max = count_max
      )
      all_paths[[combo_idx]] <- list(
        c_val = c_val,
        eta_val = eta_val,
        v_vec = v_vec,
        res = res
      )
      combo_idx <- combo_idx + 1L
    }
  }

  # Cross-sectional MAP selection
  best_path_results <- vector("list", L)

  for (l in 1:L) {
    best_score <- -Inf
    best_params_l <- list(
      lambda0 = lambda0s[l],
      c = NA,
      eta = NA,
      s = NA,
      beta = NULL,
      sigma = NA,
      score = NA
    )

    for (path in all_paths) {
      beta_l <- path$res$beta[, l]
      sigma_l <- path$res$sigmas[l]

      if (is.na(sigma_l) || sigma_l <= 0 || anyNA(beta_l)) {
        next
      }

      current_score <- compute_log_posterior(
        y = y,
        X = X,
        beta_final = beta_l,
        sigma_final = sigma_l,
        v_vec = path$v_vec,
        lambda1 = lambda1,
        lambda0_final = lambda0s[l],
        a_s = a_s,
        b_s = b_s
      )

      if (!is.finite(current_score)) {
        next
      }

      if (current_score > best_score) {
        best_score <- current_score
        best_params_l$c <- path$c_val
        best_params_l$eta <- path$eta_val
        best_params_l$s <- path$res$s_path[l]
        best_params_l$beta <- beta_l
        best_params_l$sigma <- sigma_l
        best_params_l$score <- current_score
      }
    }

    if (is.infinite(best_score)) {
      warning(sprintf("No valid score for lambda0[%d] = %.4f", l, lambda0s[l]))
    }
    best_path_results[[l]] <- best_params_l
  }

  best_path_results
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

select_lambda0_bic <- function(ssl_object, X, y) {
  rss <- apply(
    t(do.call(rbind, purrr::transpose(ssl_object)$beta)),
    2,
    function(b) sum((y - X %*% b)^2)
  )
  df <- apply(
    t(do.call(rbind, purrr::transpose(ssl_object)$beta)),
    2,
    function(b) sum(b != 0)
  )
  bic <- n * log(rss / n) + df * log(n)
  best_idx <- which.min(bic)

  c(ssl_object[[best_idx]]$intercept, ssl_object[[best_idx]]$beta)
}

# simulation function -----------------------------------------------------
# this function runs the simulation for baseline methods
baseline_data_sim_function <- function(seed, n) {
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

  baseline_fits <- list(
    "lasso" = lasso_results,
    "horseshoe" = hs_coef
  )

  if (fixed_s == TRUE) {
    # run standard discrete spike and slab
    baseline_fits$`ss, fixed s` <- lsp_fixed_ss_gibbs_sampler(
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

    # run standard spike and slab lasso
    baseline_fits$`ssl, fixed s` <- lsp_fixed_ssl_map(
      X,
      y,
      c_space = 0,
      E_space = 0,
      weights = NULL,
      s_fixed = sparsity
    ) |>
      select_lambda0_bic(X = X, y = y)
  }
  # can also evaluate LSP with random sparsity
  if (random_s == TRUE) {
    # run standard discrete spike and slab with random s
    baseline_fits$`ss, random s` <- lsp_random_ss_gibbs_sampler(
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

    # run standard spike and slab lasso
    baseline_fits$`ssl, random s` <- lsp_random_ssl_map(
      X,
      y,
      c_space = 0,
      E_space = 0,
      weights = NULL
    ) |>
      select_lambda0_bic(X = X, y = y)
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
    # run LSP for SS with fixed sparsity
    all_fits$`lsp, fixed s` <- lsp_fixed_ss_gibbs_sampler(
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

    # run LSP for SSL with fixed sparsity
    all_fits$`lsp, ssl fixed s` <- lsp_fixed_ssl_map(
      X,
      y,
      c_space = 1,
      E_space = c(0, eta_range),
      weights = weights,
      s_fixed = sparsity
    ) |>
      select_lambda0_bic(X = X, y = y)
  }

  # can also evaluate LSP with random sparsity
  if (random_s == TRUE) {
    # run LSP with random sparsity
    all_fits$`lsp, random s` <- lsp_random_ss_gibbs_sampler(
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

    # run LSP for SSL with fixed sparsity
    all_fits$`lsp, ssl random s` <- lsp_random_ssl_map(
      X,
      y,
      c_space = 1,
      E_space = c(0, eta_range),
      weights = weights
    ) |>
      select_lambda0_bic(X = X, y = y)
  }

  # return metrics
  metrics <- map(
    all_fits,
    ~ compile_model_metrics(.x, weights, alpha, beta)
  )

  return(metrics)
}
