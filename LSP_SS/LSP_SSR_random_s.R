# LLM Sparsity Prior (LSP) — Discrete Spike-and-Slab Regression with Random Sparsity
#
# Extends the fixed-sparsity LSP sampler by placing a Beta prior on the
# global sparsity parameter s, which is updated each iteration via a
# Metropolis-Hastings step on the logit scale. The concentration parameter
# eta is sampled from a discrete grid.

source("utils.R") # load helpers

# ------------------------------------------------------------------------------
# Log Prior on gamma | s, u
#
# Computes the log probability of the inclusion vector gamma under the
# factored Bernoulli prior with inclusion probabilities theta_j = s * u_j.
#
# Arguments:
#   gamma - Binary inclusion vector (length p)
#   s     - Global sparsity scalar in (0, 1)
#   u     - Weight vector (length p); theta_j = s * u_j, clamped to (0, 1)
#
# Returns:
#   Scalar log prior probability of gamma
# ------------------------------------------------------------------------------
compute_log_prior_gamma <- function(gamma, s, u) {
  theta <- pmax(1e-8, pmin(s * u, 1 - 1e-8))
  log_prob <- sum(gamma * log(theta) + (1 - gamma) * log(1 - theta))
  log_prob
}


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
#   s        - Global sparsity scalar
#   u        - Normalized weight vector (length p); theta_j = s * u_j
#   n        - Number of observations
#
# Returns:
#   Scalar unnormalized log-posterior (constants that cancel in the
#   Metropolis acceptance ratio are omitted)
# ------------------------------------------------------------------------------
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
  n_gam <- ncol(Z)

  # Flat prior on the intercept (first column of Z); slab prior on the rest
  prior_prec <- diag(c(0, rep(1 / tau, n_gam - 1)), nrow = n_gam)
  Q <- Z_gram + prior_prec
  cholQ <- chol(Q)
  log_detQ <- 2 * sum(log(diag(cholQ)))

  model_prior <- compute_log_prior_gamma(gamma, s, u)

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
lsp_random_ss_log_posterior_ratio <- function(
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

# ------------------------------------------------------------------------------
# Discrete Spike-and-Slab Gibbs Sampler with Random Sparsity
#
# Runs a Metropolis-within-Gibbs sampler for Bayesian variable selection,
# jointly inferring the sparsity level s alongside the inclusion indicators.
# When weights are supplied, the prior is modulated by LLM-derived weights
# via a concentration parameter eta sampled from a discrete grid.
#
# Arguments:
#   X                - n x p design matrix; columns are standardized internally
#                      and coefficients are returned on the original scale
#                      (intercept added internally)
#   y                - Response vector (length n)
#   weights          - Optional LLM-derived weight vector (length p); NULL disables
#   E_space          - Grid of eta values controlling weight concentration.
#                      NULL triggers automatic grid search; 0 disables weighting.
#   eta_pi_0         - Prior probability that eta = 0; the remaining mass is
#                      split uniformly over the nonzero values of E_space. Set to
#                      0 to fix eta at a single value (e.g., E_space = 2,
#                      eta_pi_0 = 0)
#   a_sigma          - Shape hyperparameter for the inverse-gamma prior on sigma^2
#   b_sigma          - Rate hyperparameter for the inverse-gamma prior on sigma^2
#   tau              - Slab variance: beta_j ~ N(0, tau * sigma^2); the
#                      intercept alpha has a flat prior
#   a_s              - First shape parameter of the Beta(a_s, b_s) prior on s
#   b_s              - Second shape parameter of the Beta(a_s, b_s) prior on s;
#                      defaults to p following Rockova & George (2018)
#   s_proposal_sigma - Standard deviation of the random-walk proposal for s
#                      on the logit scale
#   s_max            - Value of s used only to build the automatic eta grid:
#                      the largest grid value satisfies s_max * max_j u_j < 1
#   iter             - Total number of MCMC iterations
#   burn_in          - Number of initial iterations discarded as burn-in
#   thin             - Thinning interval applied after burn-in
#   prob_add         - MH proposal probability of adding a variable
#   prob_delete      - MH proposal probability of removing a variable
#                      (swap probability = 1 - prob_add - prob_delete)
#   init_weights     - If TRUE, initialize gamma using top-weighted covariates;
#                      if FALSE, initialize using marginal correlations with y
#   return_samples   - If TRUE, return all post-burn-in draws; if FALSE, return
#                      posterior means only (reduces memory for large problems)
#
# Returns:
#   A list with components:
#     beta       - Posterior draws (or mean) of the full coefficient vector
#                  (intercept first), on the original scale of X
#     gamma      - Posterior draws (or mean) of the inclusion indicators
#     invsigma_2 - Posterior draws (or mean) of the inverse noise variance
#     eta        - Posterior draws (or mean) of the concentration parameter
#     s          - Posterior draws (or mean) of the sparsity parameter
#     accs       - Metropolis acceptance indicators for gamma (or mean rate)
#     acc_s      - Metropolis acceptance indicators for s (or mean rate)
#     p_eta_0    - Full-conditional probability of eta = 0 at each draw
# ------------------------------------------------------------------------------
lsp_random_ss_gibbs_sampler <- function(
  X,
  y,
  weights = NULL,
  E_space = NULL,
  eta_pi_0 = 0.5,
  a_sigma = 1,
  b_sigma = 1,
  tau = 1,
  a_s = 1,
  b_s = NA,
  s_proposal_sigma = 1,
  s_max = 0.01,
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
  if (prob_add <= 0 || prob_delete <= 0 || prob_add + prob_delete > 1) {
    stop("prob_add and prob_delete must be positive and sum to at most 1")
  }
  if (burn_in >= iter) {
    stop("burn_in must be smaller than iter")
  }

  X_std <- standardize_X(X)
  X <- X_std$X

  # --------------------------------------------------------------------------
  # Build the eta grid and the normalized weight matrix u_mat
  # theta_j = s * u_j, where u_j = w_j^eta / mean(w^eta)
  # --------------------------------------------------------------------------

  if (is.null(weights)) {
    # No weights supplied: fix eta = 0, which sets all u_j = 1
    E_space <- 0
  } else if (is.null(E_space)) {
    # Search for the largest eta (up to 20) such that s * u_j < 1 for all j,
    # evaluated at s = s_max
    eta_max <- 0
    step_size <- 1
    theta_bound <- s_max * max(weights)^eta_max / mean(weights^eta_max)

    while (eta_max <= 20) {
      eta_max <- eta_max + step_size
      theta_bound <- s_max * max(weights)^eta_max / mean(weights^eta_max)

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

  # Build u_mat: rows index eta, columns index covariates
  u_mat <- matrix(0, nrow = K, ncol = p)

  if (length(E_space) == 1) {
    if (E_space == 0) {
      init_weights <- FALSE
      u_mat[1, ] <- rep(1, p)
    } else {
      u_mat[1, ] <-
        (weights^E_space) /
        mean(weights^E_space)
    }
  } else {
    for (k in 1:K) {
      eta_k <- E_space[k]

      u_mat[k, ] <- (weights^eta_k) / mean((weights^eta_k))
    }
  }

  # Default b_s = p; with a_s = 1 this is the Beta(1, p) prior of Rockova & George (2018)
  if (is.na(b_s)) {
    b_s <- p
  }

  # --------------------------------------------------------------------------
  # Pre-allocate storage
  # --------------------------------------------------------------------------
  n_keep <- floor((iter - burn_in) / thin)

  if (return_samples) {
    gam_store <- matrix(0, nrow = n_keep, ncol = p)
    beta_store <- matrix(0, nrow = n_keep, ncol = p + 1)
    invsigma_2_store <- numeric(n_keep)
    acc_store <- numeric(n_keep)
    acc_s_store <- numeric(n_keep)
    eta_store <- numeric(n_keep)
    s_store <- numeric(n_keep)
    pi_eta0_store <- numeric(n_keep)
  } else {
    gam_mean <- numeric(p)
    beta_mean <- numeric(p + 1)
    invsigma_2_mean <- 0
    acc_mean <- 0
    acc_s_mean <- 0
    eta_mean <- 0
    s_mean <- 0
    p_eta_0_mean <- 0
  }

  # --------------------------------------------------------------------------
  # Initialize chain
  # --------------------------------------------------------------------------

  gam_current <- rep(0, p)
  s_current <- a_s / (a_s + b_s)

  if (init_weights) {
    # initialize with the covariates receiving the highest LLM weights
    gam_current[order(weights, decreasing = TRUE)[
      1:max(2, ceiling(s_current * p))
    ]] <- 1
  } else {
    # initialize with the covariates most correlated with y
    gam_current[order(abs(cor(X, y)), decreasing = TRUE)[
      1:max(2, ceiling(s_current * p))
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

    logacc <- lsp_random_ss_log_posterior_ratio(
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
      u_mat[eta_idx_current, ],
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

    # --- Metropolis-Hastings step for s | gamma, eta ---
    # Proposal: random walk on the logit scale to respect the (0, 1) support.
    # The Jacobian of the logit transformation is included in the acceptance ratio.

    logit_s <- log(s_current / (1 - s_current))
    logit_s_new <- rnorm(1, mean = logit_s, sd = s_proposal_sigma)
    s_new <- 1 / (1 + exp(-logit_s_new))

    log_prior_s_new <- dbeta(s_new, a_s, b_s, log = TRUE)
    log_prior_s_old <- dbeta(s_current, a_s, b_s, log = TRUE)

    log_lik_s_new <- compute_log_prior_gamma(
      gam_current,
      s_new,
      u_mat[eta_idx_current, ]
    )
    log_lik_s_old <- compute_log_prior_gamma(
      gam_current,
      s_current,
      u_mat[eta_idx_current, ]
    )

    # Jacobian: |d(logit s)/ds| = 1 / [s(1-s)], so log|J| = log(s) + log(1-s)
    log_jacobian_new <- log(s_new) + log(1 - s_new)
    log_jacobian_old <- log(s_current) + log(1 - s_current)

    log_ratio_s <- (log_lik_s_new + log_prior_s_new + log_jacobian_new) -
      (log_lik_s_old + log_prior_s_old + log_jacobian_old)

    if (log(runif(1)) < log_ratio_s) {
      s_current <- s_new
      accept_s <- TRUE
    } else {
      accept_s <- FALSE
    }

    # --- Gibbs draw for eta (discrete) | gamma, s ---
    # Unnormalized log-probabilities across the eta grid. The prior is
    # zero-inflated: P(eta = 0) = eta_pi_0, with the remaining mass split
    # evenly over the other K - 1 grid values

    W <- numeric(K)
    for (k in 1:K) {
      W[k] <- compute_log_prior_gamma(gam_current, s_current, u_mat[k, ])
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
        s_store[store_i] <- s_current
        acc_store[store_i] <- acc
        acc_s_store[store_i] <- accept_s
        pi_eta0_store[store_i] <- pi_eta[1]
      } else {
        gam_mean <- gam_mean + gam_current / n_keep
        beta_mean <- beta_mean + beta_current / n_keep
        invsigma_2_mean <- invsigma_2_mean + invsigma_2_current / n_keep
        eta_mean <- eta_mean + E_space[eta_idx_current] / n_keep
        s_mean <- s_mean + s_current / n_keep
        acc_mean <- acc_mean + acc / n_keep
        acc_s_mean <- acc_s_mean + as.numeric(accept_s) / n_keep
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
      s = s_store,
      accs = acc_store,
      acc_s = acc_s_store
    )
  } else {
    list(
      beta = beta_mean,
      gamma = gam_mean,
      invsigma_2 = invsigma_2_mean,
      eta = eta_mean,
      p_eta_0 = p_eta_0_mean,
      s = s_mean,
      accs = acc_mean,
      acc_s = acc_s_mean
    )
  }
}
