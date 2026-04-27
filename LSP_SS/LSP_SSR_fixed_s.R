# LLM Sparsity Prior (with fixed prior sparsity)
# Discrete Spike and Slab for Regression with a unique constant prior inclusion probability for each covariate

# computes unnormalized log-posterior density (marginalizing out beta and sigma)
# (Z is (1, X)), does not include some constants that are cancelled in log_acceptance_rate
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
  E_space = NULL,
  sparsity = 0.05,
  a_sigma = 1,
  b_sigma = 1,
  tau = 1,
  iter = 10000,
  burn_in = 5000,
  thin = 1,
  prob_add = 1 / 3,
  prob_delete = 1 / 3,
  init_weights = TRUE,
  return_samples = TRUE
) {
  if (is.null(weights)) {
    # if no weights, then assign zero E_space
    E_space <- 0
  } else if (is.null(E_space)) {
    eta_max <- 0
    step_size <- 1
    # Calculate initial bound to ensure it starts < 1
    theta_bound <- s * max(weights)^eta_max / mean(weights^eta_max)

    while (eta_max <= 20) {
      eta_max <- eta_max + step_size # step forward
      theta_bound <- s * max(weights)^eta_max / mean(weights^eta_max)

      # if threshold is crossed, backtrack with smaller steps
      if (theta_bound >= 1) {
        # Step back to the last safe value
        eta_max <- eta_max - step_size

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
    E_space <- seq(0, eta_max, length.out = 11)
    rm(eta_max)
  }

  p <- ncol(X)
  n <- nrow(X)

  # number of values of eta in grid
  K <- length(E_space)

  # create space for theta_mat
  theta_mat <- matrix(0, nrow = K, ncol = p)

  # eta is fixed
  if (length(E_space) == 1) {
    if (E_space == 0) {
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
        (weights^E_space) /
        mean(weights^E_space)

      theta_mat[1, ] <- pmin(pmax(raw_theta, 1e-4), 1 - 1e-4)
    }
  } else {
    # eta are random
    for (k in 1:K) {
      eta_k <- E_space[k]

      raw_theta <- sparsity *
        (weights^eta_k) /
        mean((weights^eta_k))

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
  } else {
    # reduction for memory: only store means
    gam_mean <- rep(0, p)
    beta_mean <- rep(0, p + 1)
    invsigma_2_mean <- 0
    acc_mean <- 0
    eta_mean <- 0
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

  # init inv sigma^2
  invsigma_2_current <- 1 /
    (1 /
      n *
      sum((y - cbind(1, X[, which(gam_current == 1)]) %*% beta_current)^2))

  # find current index of eta in discrete uniform grid
  if (K > 1) {
    eta_idx_current <- sample(K, 1)
  } else {
    eta_idx_current <- 1
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
      theta_mat[eta_idx_current, ],
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

    # draw new eta based on gamma

    # unnormalized log probability for all K states of eta
    W <- numeric(K)
    for (k in 1:K) {
      W[k] <- sum(
        gam_current *
          log(theta_mat[k, ]) +
          (1 - gam_current) * log(1 - theta_mat[k, ])
      )
      # zero inflated prior, reduce weight for everything else:
      if (k != 1) {
        W[k] <- W[k] / (K - 1)
      }
    }
    # normalize probabilities
    # log-sum-exp trick to prevent NaN underflow
    pi_eta <- exp(W - max(W)) / sum(exp(W - max(W)))
    eta_idx_current <- sample(K, 1, prob = pi_eta)

    # store parameters
    if (i > burn_in && (i - burn_in) %% thin == 0) {
      if (return_samples) {
        store_i <- (i - burn_in) / thin # index

        gam_store[store_i, ] <- gam_current
        beta_store[store_i, ] <- beta_current
        invsigma_2_store[store_i] <- invsigma_2_current
        eta_store[store_i] <- E_space[eta_idx_current]
        acc_store[store_i] <- acc
      } else {
        gam_mean <- gam_mean + gam_current / n_keep
        beta_mean <- beta_mean + beta_current / n_keep
        invsigma_2_mean <- invsigma_2_mean + invsigma_2_current / n_keep
        eta_mean <- E_space[eta_idx_current] / n_keep
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
      "accs" = acc_store
    )
  } else {
    list(
      "beta" = beta_mean,
      "gamma" = gam_mean,
      "invsigma_2" = invsigma_2_mean,
      "eta" = eta_mean,
      "accs" = acc_mean
    )
  }
}
