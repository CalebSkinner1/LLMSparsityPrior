# LLM Sparsity Prior (LSP) — Discrete Spike-and-Slab Regression
#
# Implements a Bayesian variable selection sampler where each covariate
# has its own prior inclusion probability (theta_j), informed
# by external LLM-derived weights. The concentration parameter eta, which
# governs how strongly the weights influence the prior, is either fixed or
# sampled from a discrete grid.

source("utils.R") # load helpers

# ------------------------------------------------------------------------------
# Log-Posterior (unnormalized, marginalizing over alpha, beta and sigma^2)
#
# Arguments:
#   Z        - Selected design matrix with intercept column prepended: cbind(1, X[, gamma == 1])
#   Z_gram   - Precomputed crossprod(Z)
#   y        - Response vector (length n)
#   tau      - Slab variance: beta_j ~ N(0, tau * sigma^2); the
#              intercept alpha has a flat prior
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

  # Flat prior on the intercept (first column of Z); slab prior on the rest
  prior_prec <- diag(c(0, rep(1 / tau, n_gam - 1)), nrow = n_gam)
  Q <- Z_gram + prior_prec
  cholQ <- chol(Q)
  log_detQ <- 2 * sum(log(diag(cholQ)))

  model_prior <- sum(gamma * log(theta) + (1 - gamma) * log(1 - theta))

  y_Z <- crossprod(Z, y)
  quadratic_term <- as.numeric(t(y_Z) %*% chol2inv(cholQ) %*% y_Z)

  -(n_gam - 1) /
    2 *
    log(tau) -
    .5 * log_detQ -
    ((n - 1) / 2 + a_sigma) *
      log(.5 * (sum(y^2) - quadratic_term) + b_sigma) +
    model_prior
}

# ------------------------------------------------------------------------------
# Log-Posterior Ratio
#
# Computes log[ p(gamma_new | data) / p(gamma_old | data) ] for use in
# the Metropolis-Hastings step of the sampler.
# ------------------------------------------------------------------------------
lsp_fixed_ss_log_posterior_ratio <- function(
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
# weights = NULL or E_space = 0 recovers the standard spike-and-slab sampler.
#
# Arguments:
#   X             - n x p design matrix; columns are standardized internally and
#                   coefficients are returned on the original scale
#                   (intercept added internally)
#   y             - Response vector (length n)
#   weights       - Optional LLM-derived weight vector (length p); NULL disables
#   E_space       - Grid of eta values controlling weight concentration.
#                   NULL triggers automatic grid search; 0 disables weighting.
#   eta_pi_0      - Prior probability that eta = 0; the remaining mass is split
#                   uniformly over the nonzero values of E_space. Set to 0 to
#                   fix eta at a single value (e.g., E_space = 2, eta_pi_0 = 0)
#   sparsity      - Prior expected inclusion rate, or a length-p vector of
#                   per-covariate inclusion probabilities
#   a_sigma       - Shape hyperparameter for the inverse-gamma prior on sigma^2
#   b_sigma       - Rate hyperparameter for the inverse-gamma prior on sigma^2
#   tau           - Slab variance: beta_j ~ N(0, tau * sigma^2); the
#                   intercept alpha has a flat prior
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
#                  (intercept first), on the original scale of X
#     gamma      - Posterior draws (or mean) of the inclusion indicators
#     invsigma_2 - Posterior draws (or mean) of the inverse noise variance
#     eta        - Posterior draws (or mean) of the concentration parameter
#     accs       - Metropolis acceptance indicators (or mean acceptance rate)
#     p_eta_0    - Full-conditional probability of eta = 0 at each draw
# ------------------------------------------------------------------------------
lsp_fixed_ss_gibbs_sampler <- function(
  X,
  y,
  weights = NULL,
  E_space = NULL,
  eta_pi_0 = 0.5,
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
  # Validate inputs and standardize X
  # --------------------------------------------------------------------------
  X <- as.matrix(X)
  y <- as.numeric(y)

  if (!is.numeric(X)) {
    stop("X must be a numeric matrix")
  }
  if (length(y) != nrow(X)) {
    stop("y must have length equal to nrow(X) (", nrow(X), ")")
  }
  if (anyNA(X) || anyNA(y)) {
    stop("Missing data (NA's) detected. Eliminate missing data before calling.")
  }
  if (!is.null(weights)) {
    if (length(weights) != ncol(X)) {
      stop("weights must have length equal to ncol(X) (", ncol(X), ")")
    }
    if (any(!is.finite(weights)) || any(weights <= 0)) {
      stop("weights must be positive and finite")
    }
  }
  if (!is.null(E_space) && any(E_space < 0)) {
    stop("E_space must be nonnegative")
  }
  if (eta_pi_0 < 0 || eta_pi_0 > 1) {
    stop("eta_pi_0 must be in [0, 1]")
  }
  if (!length(sparsity) %in% c(1L, ncol(X))) {
    stop("sparsity must be length 1 or length ncol(X) (", ncol(X), ")")
  }
  if (any(sparsity <= 0) || any(sparsity >= 1)) {
    stop("all sparsity values must be strictly between 0 and 1")
  }
  if (prob_add <= 0 || prob_delete <= 0 || prob_add + prob_delete > 1) {
    stop("prob_add and prob_delete must be positive and sum to at most 1")
  }
  if (burn_in >= iter) {
    stop("burn_in must be smaller than iter")
  }

  X_std <- standardize_X(X)
  X <- X_std$X

  # --------------------------------------------------------------------------
  # Build the eta grid and corresponding prior inclusion probability matrix
  # --------------------------------------------------------------------------
  if (is.null(weights)) {
    # No weights supplied: fix eta = 0 (inclusion probabilities given by sparsity)
    E_space <- 0
  } else if (is.null(E_space)) {
    # Search for the largest eta (up to 20) such that all theta_j remain below 1
    eta_max <- 0
    step_size <- 1
    theta_bound <- max(sparsity * weights^eta_max) / mean(weights^eta_max)

    while (eta_max <= 20) {
      eta_max <- eta_max + step_size # step forward
      theta_bound <- max(sparsity * weights^eta_max) / mean(weights^eta_max)

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
    eta_max <- min(eta_max, 20)
    E_space <- seq(0, eta_max, length.out = 11)
    rm(eta_max)
  } else {
    # User-supplied grid: ensure eta = 0 is present and listed first, since
    # the zero-inflated prior mass and p_eta_0 both attach to E_space[1]
    E_space <- sort(unique(c(0, E_space)))
  }

  p <- ncol(X)
  n <- nrow(X)
  K <- length(E_space)

  # Prior inclusion probabilities: rows index eta, columns index covariates
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

      theta_mat[1, ] <- pmin(pmax(raw_theta, 1e-8), 1 - 1e-8)
    }
  } else {
    for (k in 1:K) {
      eta_k <- E_space[k]

      raw_theta <- sparsity *
        (weights^eta_k) /
        mean((weights^eta_k))

      capped_theta <- pmin(pmax(raw_theta, 1e-8), 1 - 1e-8)

      theta_mat[k, ] <- capped_theta
    }
  }

  # --------------------------------------------------------------------------
  # Pre-allocate storage
  # --------------------------------------------------------------------------

  n_keep <- floor((iter - burn_in) / thin)

  if (return_samples) {
    gam_store <- matrix(0, nrow = n_keep, ncol = p)
    beta_store <- matrix(0, nrow = n_keep, ncol = p + 1)
    invsigma_2_store <- rep(0, n_keep)
    acc_store <- rep(0, n_keep)
    eta_store <- rep(0, n_keep)
    pi_eta0_store <- numeric(n_keep)
  } else {
    gam_mean <- rep(0, p)
    beta_mean <- rep(0, p + 1)
    invsigma_2_mean <- 0
    acc_mean <- 0
    eta_mean <- 0
    p_eta_0_mean <- 0
  }

  # Initialize gamma
  gam_current <- rep(0, p)

  if (init_weights) {
    # Initialize with the covariates receiving the highest LLM weights
    gam_current[order(weights, decreasing = TRUE)[
      1:max(2, ceiling(mean(sparsity) * p))
    ]] <- 1
  } else {
    # initialize with the covariates most correlated with y
    gam_current[order(abs(cor(X, y)), decreasing = TRUE)[
      1:max(2, ceiling(mean(sparsity) * p))
    ]] <- 1
  }

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
    log_prop_ratio <- 0
    current_model_size <- sum(gam_prop)

    if (
      length(removed_gam) == 0 ||
        (length(selected_gam) > 0 && unif_gam < prob_delete)
    ) {
      # Delete a randomly chosen active variable
      chosen <- pick_one(selected_gam)
      gam_prop[chosen] <- 0
      prob_fwd <- if (current_model_size == p) 1 else prob_delete
      prob_rev <- if (current_model_size == 1) 1 else prob_add
      log_prop_ratio <- log(prob_rev) -
        log(prob_fwd) +
        log(current_model_size) -
        log(p - current_model_size + 1)
    } else if (length(selected_gam) == 0 || unif_gam < prob_delete + prob_add) {
      # Add a randomly chosen inactive variable
      chosen <- pick_one(removed_gam)
      gam_prop[chosen] <- 1
      prob_fwd <- if (current_model_size == 0) 1 else prob_add
      prob_rev <- if (current_model_size == p - 1) 1 else prob_delete
      log_prop_ratio <- log(prob_rev) -
        log(prob_fwd) +
        log(p - current_model_size) -
        log(current_model_size + 1)
    } else {
      # Swap one active and one inactive variable
      chosen1 <- pick_one(removed_gam)
      chosen2 <- pick_one(selected_gam)
      gam_prop[chosen1] <- 1
      gam_prop[chosen2] <- 0
    }

    selected_gam <- which(gam_prop == 1)
    Z_new <- cbind(1, X[, selected_gam])
    Z_new_gram <- crossprod(Z_new)

    # --- Metropolis step for gamma ---

    logacc <- lsp_fixed_ss_log_posterior_ratio(
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

    # --- Joint draw of (sigma^2, alpha, beta_gamma) | gamma, y ---
    # sigma^{-2} is drawn with (alpha, beta_gamma) integrated out, then
    # (alpha, beta_gamma) is drawn given sigma^2. This is an exact draw from
    # p(sigma^2, alpha, beta_gamma | gamma, y), so the sweep leaves the joint
    # posterior invariant after the gamma update.

    n_coef <- ncol(Z_active)
    prior_prec <- diag(c(0, rep(1 / tau, n_coef - 1)), nrow = n_coef)
    chol_Q <- chol(Z_gram_active + prior_prec)
    Zty <- crossprod(Z_active, y)
    coef_mean <- backsolve(chol_Q, forwardsolve(t(chol_Q), Zty)) # Q^{-1} Z'y
    rss_marginal <- sum(y^2) - sum(Zty * coef_mean)

    # --- Gibbs draw for sigma^{-2} | gamma, y ---

    invsigma_2_current <- rgamma(
      1,
      shape = (n - 1) / 2 + a_sigma,
      rate = b_sigma + rss_marginal / 2
    )

    # --- Gibbs draw for (alpha, beta_gamma) | sigma^2, gamma, y ---

    beta_gamma <- as.vector(
      coef_mean + backsolve(chol_Q, rnorm(n_coef)) / sqrt(invsigma_2_current)
    )

    beta_current <- numeric(p + 1)
    beta_current[c(1, active_idx + 1)] <- beta_gamma

    # --- Gibbs draw for eta (discrete) | gamma ---
    # Unnormalized log-probabilities across the eta grid. The prior is
    # zero-inflated: P(eta = 0) = eta_pi_0, with the remaining mass split
    # evenly over the other K - 1 grid values

    W <- numeric(K)
    for (k in 1:K) {
      W[k] <- sum(
        gam_current *
          log(theta_mat[k, ]) +
          (1 - gam_current) * log(1 - theta_mat[k, ])
      )
      if (K > 1) {
        W[k] <- W[k] +
          if (k == 1) log(eta_pi_0) else log(1 - eta_pi_0) - log(K - 1)
      }
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
        pi_eta0_store[store_i] <- pi_eta[1]
      } else {
        gam_mean <- gam_mean + gam_current / n_keep
        beta_mean <- beta_mean + beta_current / n_keep
        invsigma_2_mean <- invsigma_2_mean + invsigma_2_current / n_keep
        eta_mean <- eta_mean + E_space[eta_idx_current] / n_keep
        acc_mean <- acc_mean + acc / n_keep
        p_eta_0_mean <- p_eta_0_mean + pi_eta[1] / n_keep
      }
    }
  }

  # Map coefficients back to the original scale of X
  if (return_samples) {
    beta_store <- unstandardize_beta(beta_store, X_std$center, X_std$scale)
  } else {
    beta_mean <- unstandardize_beta(beta_mean, X_std$center, X_std$scale)
  }

  if (return_samples) {
    list(
      beta = beta_store,
      gamma = gam_store,
      invsigma_2 = invsigma_2_store,
      eta = eta_store,
      p_eta_0 = pi_eta0_store,
      accs = acc_store
    )
  } else {
    list(
      beta = beta_mean,
      gamma = gam_mean,
      invsigma_2 = invsigma_2_mean,
      eta = eta_mean,
      p_eta_0 = p_eta_0_mean,
      accs = acc_mean
    )
  }
}
