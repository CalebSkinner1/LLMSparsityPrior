# LLM Sparsity Prior (LSP) — Discrete Spike-and-Slab Regression
#
# Implements a Bayesian variable selection sampler where each covariate
# has a unique prior inclusion probability (theta_j), optionally informed
# by external LLM-derived weights. Supports a fixed or grid-searched
# concentration parameter eta that governs how strongly the weights
# influence the prior.

# ------------------------------------------------------------------------------
# Log-Posterior (unnormalized, marginalizing over beta and sigma^2)
#
# Arguments:
#   Z        - Selected design matrix with intercept column prepended: cbind(1, X[, gamma == 1])
#   Z_gram   - Precomputed crossprod(Z)
#   y        - Response vector (length n)
#   tau      - Slab variance for the g-prior
#   gamma    - Binary inclusion vector (length p)
#   a_sigma  - Shape hyperparameter for the inverse-gamma prior on sigma^2
#   b_sigma  - Rate hyperparameter for the inverse-gamma prior on sigma^2
#   theta    - Prior inclusion probability vector (length p)
#   n        - Number of observations
#
# Returns:
#   Scalar unnormalized log-posterior (constants that cancel in the
#   Metropolis acceptance ratio are omitted)
# ------------------------------------------------------------------------------
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
  n_gam <- ncol(Z)

  Q <- Z_gram + diag(n_gam) / tau
  cholQ <- chol(Q)
  log_detQ <- 2 * sum(log(diag(cholQ)))

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

# ------------------------------------------------------------------------------
# Metropolis Log-Acceptance Rate
#
# Computes log[ p(gamma_new | data) / p(gamma_old | data) ] for use in
# the Metropolis-Hastings step of the sampler.
# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
# Discrete Spike-and-Slab Gibbs Sampler
#
# Runs a Metropolis-within-Gibbs sampler for Bayesian variable selection.
# When weights are supplied, the prior inclusion probabilities are modulated
# by an LLM-derived weight vector via a concentration parameter eta. Setting
# E_space = 0 (no weights) recovers the standard spike-and-slab sampler.
#
# Arguments:
#   X             - n x p design matrix (uncentered; intercept added internally)
#   y             - Response vector (length n)
#   weights       - Optional LLM-derived weight vector (length p); NULL disables
#   E_space       - Grid of eta values controlling weight concentration.
#                   NULL triggers automatic grid search; 0 disables weighting.
#   sparsity      - Prior expected inclusion rate, or a length-p vector of
#                   per-covariate inclusion probabilities
#   a_sigma       - Shape hyperparameter for the inverse-gamma prior on sigma^2
#   b_sigma       - Rate hyperparameter for the inverse-gamma prior on sigma^2
#   tau           - Slab variance for the g-prior on regression coefficients
#   iter          - Total number of MCMC iterations
#   burn_in       - Number of initial iterations discarded as burn-in
#   thin          - Thinning interval applied after burn-in
#   prob_add      - MH proposal probability of adding a variable
#   prob_delete   - MH proposal probability of removing a variable
#                   (swap probability = 1 - prob_add - prob_delete)
#   init_weights  - If TRUE, initialize gamma using top-weighted covariates;
#                   if FALSE, initialize using marginal correlations with y
#   return_samples - If TRUE, return all post-burn-in draws; if FALSE, return
#                   posterior means only (reduces memory for large problems)
#
# Returns:
#   A list with components:
#     beta       - Posterior draws (or mean) of the full coefficient vector
#     gamma      - Posterior draws (or mean) of the inclusion indicators
#     invsigma_2 - Posterior draws (or mean) of the inverse noise variance
#     eta        - Posterior draws (or mean) of the concentration parameter
#     accs       - Metropolis acceptance indicators (or mean acceptance rate)
# ------------------------------------------------------------------------------
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
  # --------------------------------------------------------------------------
  # Build the eta grid and corresponding prior inclusion probability matrix
  # --------------------------------------------------------------------------
  if (is.null(weights)) {
    # No weights supplied: fix eta = 0 (uniform inclusion probability)
    E_space <- 0
  } else if (is.null(E_space)) {
    # Search for the largest eta such that all theta_j remain below 1
    eta_max <- 0
    step_size <- 1
    theta_bound <- sparsity * max(weights)^eta_max / mean(weights^eta_max)

    while (eta_max <= 20) {
      eta_max <- eta_max + step_size # step forward
      theta_bound <- sparsity * max(weights)^eta_max / mean(weights^eta_max)

      if (theta_bound >= 1) {
        eta_max <- eta_max - step_size

        if (step_size == 1) {
          step_size <- 0.1
        } else if (step_size == 0.1) {
          step_size <- 0.01
        } else {
          break
        }
      }
    }
    E_space <- seq(0, eta_max, length.out = 11)
    rm(eta_max)
  }

  p <- ncol(X)
  n <- nrow(X)
  K <- length(E_space)

  # Compute the p(gamma_j = 1) matrix: rows index eta, columns index covariates
  theta_mat <- matrix(0, nrow = K, ncol = p)

  if (length(E_space) == 1) {
    if (E_space == 0) {
      # No weight modulation: use sparsity directly as inclusion probabilities
      init_weights <- FALSE
      if (length(sparsity) == 1) {
        theta_mat[1, ] <- rep(sparsity, p)
      } else {
        theta_mat[1, ] <- sparsity
      }
    } else {
      raw_theta <- sparsity *
        (weights^E_space) /
        mean(weights^E_space)

      theta_mat[1, ] <- pmin(pmax(raw_theta, 1e-4), 1 - 1e-4)
    }
  } else {
    for (k in 1:K) {
      eta_k <- E_space[k]

      raw_theta <- sparsity *
        (weights^eta_k) /
        mean((weights^eta_k))

      capped_theta <- pmin(pmax(raw_theta, 1e-4), 1 - 1e-4)

      theta_mat[k, ] <- capped_theta
    }
  }

  # --------------------------------------------------------------------------
  # Pre-allocate storage
  # --------------------------------------------------------------------------

  n_keep <- ceiling((iter - burn_in) / thin)

  if (return_samples) {
    gam_store <- matrix(0, nrow = n_keep, ncol = p)
    beta_store <- matrix(0, nrow = n_keep, ncol = p + 1)
    invsigma_2_store <- rep(0, n_keep)
    acc_store <- rep(0, n_keep)
    eta_store <- rep(0, n_keep)
  } else {
    gam_mean <- rep(0, p)
    beta_mean <- rep(0, p + 1)
    invsigma_2_mean <- 0
    acc_mean <- 0
    eta_mean <- 0
  }

  # generate initial values of gamma (in a smart way)
  gam_current <- rep(0, p)

  if (init_weights) {
    # Initialize with the covariates receiving the highest LLM weights
    gam_current[order(weights, decreasing = TRUE)[
      1:max(2, ceiling(sparsity * p))
    ]] <- 1
  } else {
    # initialize with the covariates most correlated with y
    gam_current[order(abs(cor(X, y)), decreasing = TRUE)[
      1:max(2, ceiling(sparsity * p))
    ]] <- 1
  }

  # Initialize beta via ridge regression on the selected covariates
  fit <- glmnet::glmnet(
    X[, which(gam_current == 1)],
    y,
    alpha = 0,
    lambda = 1
  )
  beta_current <- as.vector(coef(fit))

  # Initialize inverse noise variance from current residual variance
  invsigma_2_current <- 1 /
    (1 /
      n *
      sum((y - cbind(1, X[, which(gam_current == 1)]) %*% beta_current)^2))

  eta_idx_current <- if (K > 1) sample(K, 1) else 1

  # --------------------------------------------------------------------------
  # Main MCMC Loop
  # --------------------------------------------------------------------------
  for (i in 1:iter) {
    # --- Propose a new gamma via add / delete / swap (ADS) ---

    gam_prop <- gam_current
    selected_gam <- which(gam_prop == 1)
    removed_gam <- which(gam_prop == 0)

    Z_old <- cbind(1, X[, selected_gam])
    Z_old_gram <- crossprod(Z_old)

    unif_gam <- runif(1)
    if (length(selected_gam) == 0) {
      unif_gam <- 0.5
    } # force an add when model is empty

    log_prop_ratio <- 0
    current_model_size <- sum(gam_prop)

    if (unif_gam < prob_delete || length(removed_gam) == 0) {
      # Delete a randomly chosen active variable
      chosen <- sample(selected_gam)[1]
      gam_prop[chosen] <- 0
      log_prop_ratio <- log(prob_add) -
        log(prob_delete) +
        log(current_model_size) -
        log(p - current_model_size + 1)
    } else if (unif_gam < prob_delete + prob_add) {
      # Add a randomly chosen inactive variable
      chosen <- sample(removed_gam)[1]
      gam_prop[chosen] <- 1
      log_prop_ratio <- log(prob_delete) -
        log(prob_add) +
        log(p - current_model_size) -
        log(current_model_size + 1)
    } else {
      # Swap one active and one inactive variable
      chosen1 <- sample(removed_gam)[1]
      chosen2 <- sample(selected_gam)[1]
      gam_prop[chosen1] <- 1
      gam_prop[chosen2] <- 0
    }

    selected_gam <- which(gam_prop == 1)
    Z_new <- cbind(1, X[, selected_gam])
    Z_new_gram <- crossprod(Z_new)

    # --- Metropolis step for gamma ---

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
      gam_current <- gam_prop
      Z_active <- Z_new
      Z_gram_active <- Z_new_gram
      active_idx <- selected_gam
      acc <- 1
    } else {
      Z_active <- Z_old
      Z_gram_active <- Z_old_gram
      active_idx <- which(gam_current == 1)
      acc <- 0
    }

    # --- Gibbs draw for beta | gamma, sigma^2 ---

    chol_mat <- chol(
      invsigma_2_current *
        Z_gram_active +
        (invsigma_2_current / tau) * diag(ncol(Z_gram_active))
    )
    invQ <- chol2inv(chol_mat)
    l <- invsigma_2_current * crossprod(Z_active, y)
    beta_gamma <- MASS::mvrnorm(1, invQ %*% l, invQ)

    beta_current <- numeric(p + 1)
    beta_current[c(1, active_idx + 1)] <- beta_gamma

    # --- Gibbs draw for sigma^{-2} | gamma, beta ---

    invsigma_2_current <- rgamma(
      1,
      shape = n / 2 + a_sigma,
      rate = 0.5 * sum((y - Z_active %*% beta_gamma)^2) + b_sigma
    )

    # --- Gibbs draw for eta (discrete) | gamma ---
    # Unnormalized log-probabilities across the eta grid; eta = 0 receives
    # a zero-inflated weight (down-weighted by 1/(K-1) for all other states)

    W <- numeric(K)
    for (k in 1:K) {
      W[k] <- sum(
        gam_current *
          log(theta_mat[k, ]) +
          (1 - gam_current) * log(1 - theta_mat[k, ])
      )
      if (k != 1) W[k] <- W[k] / (K - 1)
    }
    pi_eta <- exp(W - max(W)) / sum(exp(W - max(W))) # log-sum-exp stabilization
    eta_idx_current <- sample(K, 1, prob = pi_eta)

    # --- Store post-burn-in draws ---

    if (i > burn_in && (i - burn_in) %% thin == 0) {
      if (return_samples) {
        store_i <- (i - burn_in) / thin
        gam_store[store_i, ] <- gam_current
        beta_store[store_i, ] <- beta_current
        invsigma_2_store[store_i] <- invsigma_2_current
        eta_store[store_i] <- E_space[eta_idx_current]
        acc_store[store_i] <- acc
      } else {
        gam_mean <- gam_mean + gam_current / n_keep
        beta_mean <- beta_mean + beta_current / n_keep
        invsigma_2_mean <- invsigma_2_mean + invsigma_2_current / n_keep
        eta_mean <- eta_mean + E_space[eta_idx_current] / n_keep
        acc_mean <- acc_mean + acc / n_keep
      }
    }
  }

  if (return_samples) {
    list(
      beta = beta_store,
      gamma = gam_store,
      invsigma_2 = invsigma_2_store,
      eta = eta_store,
      accs = acc_store
    )
  } else {
    list(
      beta = beta_mean,
      gamma = gam_mean,
      invsigma_2 = invsigma_2_mean,
      eta = eta_mean,
      accs = acc_mean
    )
  }
}
