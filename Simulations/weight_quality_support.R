# Simulation Helper Functions
#
# Utility functions supporting the simulation drivers (weight_quality_sims.R,
# block_corr_sims.R, and eta_sensitivity_sims.R). Provides:
#   - Weight agreement metrics (L1, L2, pairwise, ROC-AUC)
#   - Synthetic weight generation from a target L1 agreement level
#   - Model fit metric compilation
#   - Three simulation entry points:
#       baseline_data_sim_function  -- fits weight-free baseline models
#       sim_function                -- fits LSP models given a weight vector
#       eta_sensitivity_function    -- fits LSP models at a fixed eta value
#
# Note: the simulation functions consume several global variables defined in
# the calling driver script (e.g., p, s, n, sparsity, tau, iter, burn_in,
# effect_size, y_sd, Xvar, cov_mat, eta_range, fixed_s, random_s, a_sigma,
# b_sigma). They must be sourced after those globals are set.

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library("dplyr")
  library("tidyr")
  library("purrr")
  library("stringr")
  library("readr")
  library("glmnet")
  library("parallel")
  library("furrr")
})

source("utils.R") # shared helpers, including llm_lasso_simp
source("LSP_SS/LSP_SSR_fixed_s.R")
source("LSP_SS/LSP_SSR_random_s.R")
source("LSP_SSL/LSP_SSLR.R")


# ------------------------------------------------------------------------------
# Weight Agreement Metrics
#
# Each metric compares a binary ground-truth inclusion vector (true_gamma)
# to a continuous weight vector, returning a scalar where higher values
# indicate better agreement.
# ------------------------------------------------------------------------------

# L1 agreement: 1 - mean absolute deviation after min-max scaling
l1_weight_agreement <- function(true_gamma, weights) {
  scaled_weights <- (weights - min(weights)) / (max(weights) - min(weights))
  1 - mean(abs(true_gamma - scaled_weights))
}

# L2 agreement: 1 - mean squared deviation after min-max scaling
l2_weight_agreement <- function(true_gamma, weights) {
  scaled_weights <- (weights - min(weights)) / (max(weights) - min(weights))
  1 - mean((true_gamma - scaled_weights)^2)
}

# Pairwise agreement: fraction of covariate pairs ranked
# consistently between true_gamma and weights
pairwise_weight_agreement <- function(true_gamma, weights) {
  if (length(true_gamma) != length(weights)) {
    stop("true_gamma and weights must have the same number of elements")
  }

  true_gamma_mat <- ifelse(outer(true_gamma, true_gamma, FUN = "-") > 0, 1, 0)
  weights_mat <- ifelse(outer(weights, weights, FUN = "-") > 0, 1, 0)

  total_disagreement <- sum(abs(true_gamma_mat - weights_mat))
  total_pairs <- length(true_gamma) * (length(true_gamma) - 1) / 2

  (total_pairs - total_disagreement) / total_pairs
}

# ROC-AUC agreement: area under the ROC curve treating weights as a classifier
# for true_gamma (requires pROC)
ROC_weight_agreement <- function(true_gamma, weights) {
  scaled_weights <- (weights - min(weights)) / (max(weights) - min(weights))
  ROC <- pROC::roc(
    true_gamma,
    scaled_weights,
    levels = c("0", "1"),
    direction = "<"
  )
  pROC::auc(ROC)[1]
}

# ------------------------------------------------------------------------------
# Block Correlation Structure
#
# Constructs a p x p block-diagonal correlation matrix with equal-sized blocks.
# Within each block, covariates share a common (exchangeable) correlation rho;
# covariates in different blocks are uncorrelated.
#
# Arguments:
#   p          - Number of covariates
#   block_size - Number of covariates per block (must divide p)
#   rho        - Within-block correlation
#
# Returns:
#   A p x p correlation matrix
# ------------------------------------------------------------------------------
block_cor_mat <- function(p, block_size, rho) {
  if (p %% block_size != 0) {
    stop("p must be divisible by block_size")
  }
  if (rho <= -1 / (block_size - 1) || rho >= 1) {
    stop("rho must lie in (-1 / (block_size - 1), 1) for positive definiteness")
  }

  block <- matrix(rho, nrow = block_size, ncol = block_size)
  diag(block) <- 1

  kronecker(diag(p / block_size), block)
}


# ------------------------------------------------------------------------------
# Synthetic Weight Generation
#
# Constructs an integer-valued weight vector with a specified L1 agreement
# (phi) relative to true_gamma. Weights for signal covariates (gamma = 1) are
# drawn from the high end of a discrete scale [1, categories]; weights for
# noise covariates (gamma = 0) are drawn from the low end. The proportions
# across rank categories follow a geometric series whose ratio is solved
# numerically to achieve the target phi.
#
# Arguments:
#   phi        - Target L1 weight agreement in [0.5, 1]; phi = 0.5 gives
#                uninformative weights and phi = 1 gives perfect weights
#   true_gamma - Binary ground-truth inclusion vector
#   categories - Number of discrete weight levels (default 5)
#
# Returns:
#   Integer weight vector (length p)
# ------------------------------------------------------------------------------
generate_weights <- function(phi, true_gamma, categories = 5) {
  if (phi < 0.5 || phi > 1) {
    stop("phi must be in [0.5, 1]")
  }

  if (phi == 1) {
    weight_prop <- c(1, rep(0, categories - 1))
  } else {
    mu <- (categories - 1) * (1 - phi)

    # Solve for the geometric ratio r such that the induced mean equals mu
    f <- function(r) {
      k <- 1:(categories - 1)
      sum((k - mu) * r^k) - mu
    }
    r <- uniroot(f, interval = c(0, 1))$root

    k <- 1:(categories - 1)
    c_scale <- mu / sum(k * r^k)
    weight_prop <- c_scale * (r^(0:(categories - 1)))
  }

  # Distribute counts across rank categories, correcting rounding remainders
  assign_counts <- function(group_size) {
    dist <- weight_prop * group_size
    counts <- floor(dist)
    leftover <- group_size - sum(counts)
    if (leftover > 0) {
      extra <- order(dist - counts, decreasing = TRUE)[1:leftover]
      counts[extra] <- counts[extra] + 1
    }
    counts
  }

  # Noise covariates receive low ranks; signal covariates receive high ranks
  counts_0 <- assign_counts(sum(true_gamma == 0))
  counts_1 <- assign_counts(sum(true_gamma == 1))

  weights_0 <- rep(seq_along(counts_0), times = round(counts_0))
  weights_1 <- rep(
    (categories + 1) - seq_along(counts_1),
    times = round(counts_1)
  )

  weights <- numeric(length(true_gamma))
  weights[true_gamma == 0] <- weights_0
  weights[true_gamma == 1] <- weights_1

  weights
}


# ------------------------------------------------------------------------------
# Model Fit Metric Compilation
#
# Extracts selection indicators and coefficient estimates from a fitted model
# object and computes classification metrics (F1, FP, FN) and L1 coefficient
# error. Handles both MCMC output (list with $gamma and $beta) and MAP/Lasso
# output (numeric coefficient vector including intercept).
#
# Arguments:
#   model_object - MCMC list (with $gamma, $beta) or numeric coefficient vector
#   weights      - Weight vector used in fitting (length p; for stratified summaries)
#   alpha        - True intercept scalar
#   beta         - True coefficient vector (length p)
#
# Returns:
#   A one-row tibble with columns for per-weight-level mean selection rates
#   (group_type) plus f1, l1, fp, and fn
# ------------------------------------------------------------------------------
compile_model_metrics <- function(model_object, weights, alpha, beta) {
  signal <- beta != 0

  if (is.numeric(model_object)) {
    gamma_predict <- as.vector(model_object != 0)[-1]
    l1 <- sum(abs(as.vector(model_object) - c(alpha, beta)))
    p_eta_0 <- NA_real_
    eta_post <- NA_real_
  } else {
    gamma_predict <- model_object$gamma
    l1 <- sum(abs(model_object$beta - c(alpha, beta)))
    p_eta_0 <- model_object$p_eta_0 %||% NA_real_
    eta_post <- model_object$eta %||% NA_real_
  }

  fp <- length(which(gamma_predict[!signal] > 0.5))
  fn <- length(which(gamma_predict[signal] < 0.5))
  tp <- sum(signal) - fn
  precision <- if_else(tp == 0, 0, tp / (tp + fp))
  recall <- if_else(tp == 0, 0, tp / (tp + fn))

  tibble(gamma_predict, weights, signal) |>
    group_by(weights, signal) |>
    summarize(coef_mean = mean(gamma_predict), .groups = "keep") |>
    arrange(signal) |>
    mutate(group_type = str_c("w", weights, "s", as.numeric(signal))) |>
    ungroup() |>
    select(group_type, coef_mean) |>
    pivot_wider(names_from = group_type, values_from = coef_mean) |>
    mutate(
      f1 = if_else(
        precision + recall == 0,
        0,
        2 * (precision * recall) / (precision + recall)
      ),
      l1 = l1,
      fp = fp,
      fn = fn,
      p_eta_0 = p_eta_0,
      eta_post = eta_post
    )
}

# ------------------------------------------------------------------------------
# Simulation Functions
#
# All three functions below read the following globals from the calling driver
# script (weight_quality_sims.R, block_corr_sims.R, or eta_sensitivity_sims.R):
#   p, s, n, effect_size, y_sd, Xvar, cov_mat, sparsity, a_sigma, b_sigma,
#   tau, iter, burn_in, eta_range, fixed_s, random_s
# ------------------------------------------------------------------------------

# Generate data and fit weight-free baseline models (Lasso, horseshoe, and
# optionally standard SS / SSL without LLM weights). Returns a list containing
# the generated dataset and fitted baseline objects, to be passed to sim_function.
baseline_data_sim_function <- function(seed, n, randomize_beta = FALSE) {
  set.seed(seed)

  X <- MASS::mvrnorm(n, mu = rep(0, p), Xvar * cov_mat)

  # permute beta ordering
  perm <- if (randomize_beta) sample.int(p) else seq_len(p)
  beta <- c(rep(0, p - s), rep(effect_size, s))[perm]
  alpha <- effect_size
  y <- X %*% beta + alpha + rnorm(n, 0, sd = y_sd)

  lasso_cv <- glmnet::cv.glmnet(X, y, alpha = 1)
  lasso_results <- as.vector(coef(lasso_cv, s = "lambda.min"))

  hs_fit <- Mhorseshoe::approx_horseshoe(
    y = y,
    X = cbind(1, X),
    burn = 10000,
    iter = 5000
  )
  # A covariate is selected when its 95% credible interval excludes
  # [-1e-4, 1e-4]; the first element (intercept) is not a covariate
  hs_selected <- hs_fit$LeftCI > 1e-4 | hs_fit$RightCI < -1e-4
  hs_coef <- list(beta = hs_fit$BetaHat, gamma = as.numeric(hs_selected[-1]))
  rm(hs_fit)
  gc()

  baseline_fits <- list(lasso = lasso_results, horseshoe = hs_coef)

  if (fixed_s) {
    baseline_fits$`ss, fixed s` <- lsp_fixed_ss_gibbs_sampler(
      X,
      y,
      E_space = 0,
      sparsity = sparsity,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau,
      iter = iter,
      burn_in = burn_in,
      init_weights = FALSE,
      return_samples = FALSE
    )
    baseline_fits$`ssl, fixed s` <- lsp_ssl_map(
      X,
      y,
      E_space = 0,
      weights = NULL,
      penalty = "separable",
      variance = "fixed",
      sparsity = sparsity
    ) |>
      select_lambda0_bic(X = X, y = y)
  }

  if (random_s) {
    baseline_fits$`ss, random s` <- lsp_random_ss_gibbs_sampler(
      X,
      y,
      E_space = 0,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau,
      iter = iter,
      burn_in = burn_in,
      init_weights = FALSE,
      return_samples = FALSE
    )
    baseline_fits$`ssl, random s` <- lsp_ssl_map(
      X,
      y,
      E_space = 0,
      weights = NULL,
      penalty = "adaptive",
      variance = "fixed"
    ) |>
      select_lambda0_bic(X = X, y = y)
  }

  list(
    data = list(X = X, y = y, alpha = alpha, beta = beta, perm = perm),
    baselines = baseline_fits
  )
}

# Fit all LSP models (and LLM-Lasso) for a given weight vector, reusing the
# data and baseline fits produced by baseline_data_sim_function. Returns a
# named list of metric tibbles, one per method.
sim_function <- function(baseline_fits, weights) {
  X <- baseline_fits$data$X
  y <- baseline_fits$data$y
  beta <- baseline_fits$data$beta
  alpha <- baseline_fits$data$alpha

  weights <- weights[baseline_fits$data$perm]

  all_fits <- baseline_fits$baselines

  all_fits$`llm-lasso` <- llm_lasso_simp(
    X_train = X,
    y_train = y,
    weights = weights,
    elastic_net = 1,
    regression = TRUE
  )$coef

  if (fixed_s) {
    all_fits$`lsp, fixed s` <- lsp_fixed_ss_gibbs_sampler(
      X,
      y,
      weights,
      sparsity = sparsity,
      E_space = eta_range,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau,
      iter = iter,
      burn_in = burn_in,
      init_weights = TRUE,
      return_samples = FALSE
    )
    all_fits$`lsp, ssl fixed s` <- lsp_ssl_map(
      X,
      y,
      E_space = eta_range,
      weights = weights,
      penalty = "separable",
      variance = "fixed",
      sparsity = sparsity
    ) |>
      select_lambda0_bic(X = X, y = y)
  }

  if (random_s) {
    all_fits$`lsp, random s` <- lsp_random_ss_gibbs_sampler(
      X,
      y,
      weights,
      E_space = eta_range,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau,
      iter = iter,
      burn_in = burn_in,
      init_weights = TRUE,
      return_samples = FALSE
    )
    all_fits$`lsp, ssl random s` <- lsp_ssl_map(
      X,
      y,
      E_space = eta_range,
      weights = weights,
      penalty = "adaptive",
      variance = "fixed"
    ) |>
      select_lambda0_bic(X = X, y = y)
  }

  map(all_fits, ~ compile_model_metrics(.x, weights, alpha, beta))
}

# Fit LSP models at a user-specified fixed eta value (set_eta) for a single
# simulation replicate. eta is held fixed by placing zero prior mass on
# eta = 0 (eta_pi_0 = 0). Used to assess sensitivity to the choice of eta.
eta_sensitivity_function <- function(
  seed,
  weights,
  set_eta,
  randomize_beta = FALSE
) {
  set.seed(seed)

  X <- MASS::mvrnorm(n, mu = rep(0, p), Xvar * cov_mat)

  # permute beta ordering, and align weights to the same ordering
  perm <- if (randomize_beta) sample.int(p) else seq_len(p)
  beta <- c(rep(0, p - s), rep(effect_size, s))[perm]
  weights <- weights[perm]
  alpha <- effect_size
  y <- X %*% beta + alpha + rnorm(n, 0, sd = y_sd)

  all_fits <- list()

  if (fixed_s) {
    all_fits$`lsp, fixed s` <- lsp_fixed_ss_gibbs_sampler(
      X,
      y,
      weights,
      sparsity = sparsity,
      E_space = set_eta,
      eta_pi_0 = 0,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau,
      iter = iter,
      burn_in = burn_in,
      init_weights = TRUE,
      return_samples = FALSE
    )
    all_fits$`lsp, ssl fixed s` <- lsp_ssl_map(
      X,
      y,
      E_space = set_eta,
      eta_pi_0 = 0,
      weights = weights,
      penalty = "separable",
      variance = "fixed",
      sparsity = sparsity
    ) |>
      select_lambda0_bic(X = X, y = y)
  }

  if (random_s) {
    all_fits$`lsp, random s` <- lsp_random_ss_gibbs_sampler(
      X,
      y,
      weights,
      E_space = set_eta,
      eta_pi_0 = 0,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau,
      iter = iter,
      burn_in = burn_in,
      init_weights = TRUE,
      return_samples = FALSE
    )
    all_fits$`lsp, ssl random s` <- lsp_ssl_map(
      X,
      y,
      E_space = set_eta,
      eta_pi_0 = 0,
      weights = weights,
      penalty = "adaptive",
      variance = "fixed"
    ) |>
      select_lambda0_bic(X = X, y = y)
  }

  map(all_fits, ~ compile_model_metrics(.x, weights, alpha, beta))
}
