# Simulation Helper Functions
#
# Utility functions supporting weight_quality_sims.R. Provides:
#   - Weight agreement metrics (L1, L2, pairwise, ROC-AUC)
#   - Synthetic weight generation from a target L1 agreement level
#   - Model fit metric compilation
#   - LLM-Lasso implementation (adapted from Zhang et al.; see attribution below)
#   - BIC-based lambda selection for SSL outputs
#   - Three simulation entry points:
#       baseline_data_sim_function  -- fits weight-free baseline models
#       sim_function                -- fits LSP models given a weight vector
#       eta_sensitivity_function    -- fits LSP models at a fixed eta value
#
# Note: the simulation functions consume several global variables defined in
# weight_quality_sims.R (e.g., p, s, n, sparsity, tau, iter, burn_in,
# effect_size, y_sd, Xvar, cov_mat, eta_range, fixed_s, random_s, a_sigma,
# b_sigma). They must be sourced after those globals are set.

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library("tidyverse")
  library("hypergeo")
  library("glmnet")
  library("parallel")
  library("furrr")
})

source("LSP_SS/LSP_SSR_fixed_s.R")
source("LSP_SS/LSP_SSR_random_s.R")
source("LSP_SSL/LSP_SSLR.R")


# ------------------------------------------------------------------------------
# Weight Agreement Metrics
#
# Each metric compares a binary ground-truth inclusion vector (true_gamma)
# to a continuous weight vector, returning a scalar in [0, 1] where higher
# values indicate better agreement.
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
#   phi        - Target L1 weight agreement in (0, 1]; phi = 1 gives perfect weights
#   true_gamma - Binary ground-truth inclusion vector
#   categories - Number of discrete weight levels (default 5)
#
# Returns:
#   Integer weight vector (length p)
# ------------------------------------------------------------------------------
generate_weights <- function(phi, true_gamma, categories = 5) {
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
    # Lasso / MAP: coefficient vector includes intercept at position 1
    gamma_predict <- as.vector(model_object != 0)[-1]
    l1 <- sum(abs(as.vector(model_object) - c(alpha, beta)))
  } else {
    # MCMC: posterior means of gamma and beta stored in list
    gamma_predict <- model_object$gamma
    l1 <- sum(abs(model_object$beta - c(alpha, beta)))
  }

  fp <- length(which(gamma_predict[!signal] > 0.5))
  fn <- length(which(gamma_predict[signal] < 0.5))
  tp <- sum(signal) - fn
  precision <- if_else(tp == 0, 0, tp / (tp + fp))
  recall <- if_else(tp == 0, 0, tp / (tp + fn))

  tibble(gamma_predict, weights, signal) %>%
    group_by(weights, signal) %>%
    summarize(coef_mean = mean(gamma_predict), .groups = "keep") %>%
    arrange(signal) %>%
    mutate(group_type = str_c("w", weights, "s", as.numeric(signal))) %>%
    ungroup() %>%
    select(group_type, coef_mean) %>%
    pivot_wider(names_from = group_type, values_from = coef_mean) %>%
    mutate(
      f1 = if_else(
        precision + recall == 0,
        0,
        2 * (precision * recall) / (precision + recall)
      ),
      l1 = l1,
      fp = fp,
      fn = fn
    )
}

# ------------------------------------------------------------------------------
# LLM-Lasso
#
# Adapted from Zhang et al.: https://github.com/pilancilab/LLM-Lasso
# Lightly edited for compatibility with this simulation framework.
# ------------------------------------------------------------------------------

# Scale X_new using the center and standard deviation computed from X_train.
# Columns with zero or non-finite variance are left unscaled.
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

# Align a (optionally named) weight vector to X's column order and validate
# that all entries are strictly positive and finite.
.align_and_check_weights <- function(weights, X) {
  if (!is.null(names(weights))) {
    missing_cols <- setdiff(colnames(X), names(weights))
    if (length(missing_cols) > 0) {
      stop(
        "weights are named but missing entries for features: ",
        paste(missing_cols, collapse = ", ")
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

# Compute the area between the candidate CV error curve and the baseline
# (uniform penalty) CV error curve, interpolated to a common sparsity grid.
# A larger value indicates that the candidate penalty factor outperforms the
# baseline across the regularization path.
cve <- function(cvm, non_zero, ref_cvm, ref_non_zero) {
  df1 <- tibble(ref_non_zero, ref_cvm) %>%
    group_by(ref_non_zero) %>%
    summarise(ref_cvm = min(ref_cvm), .groups = "drop") %>%
    arrange(ref_non_zero)
  df2 <- tibble(non_zero, cvm) %>%
    group_by(non_zero) %>%
    summarise(cvm = min(cvm), .groups = "drop") %>%
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

# Fit LLM-Lasso by selecting the penalty factor exponent (1/w^k, k = 0,...,
# max_imp_pow) that maximizes the area between its CV error curve and the
# unweighted baseline, then re-fitting at the chosen penalty.
#
# Arguments:
#   X_train          - Training predictor matrix
#   y_train          - Training response vector or factor
#   weights          - LLM-derived weight vector (length = ncol(X_train))
#   folds_cv         - Number of cross-validation folds
#   elastic_net      - Elastic net mixing parameter (1 = lasso, 0 = ridge)
#   max_imp_pow      - Maximum exponent for the penalty factor search
#   lambda_min_ratio - Ratio of smallest to largest lambda in the path
#   regression       - If TRUE, fit Gaussian regression; else classification
#   multinomial      - If TRUE, use multinomial family (overrides regression)
#   type_measure     - CV loss metric; defaults to "mse" (regression) or
#                      "class" (classification)
#   use_lambda_1se   - If TRUE, use lambda.1se; otherwise use lambda.min
#
# Returns:
#   A list with components: model (chosen penalty name), algo, method
#   (selected lambda), coef (coefficient vector or per-class list),
#   n_features (number of non-zero predictors)
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
  glm_family <- if (multinomial) {
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
  if (glm_family %in% c("binomial", "multinomial")) {
    y_train <- if (is.factor(y_train)) y_train else factor(y_train)
  }

  X_train_sc <- .scale_like_train(X_train)$X_train
  w <- .align_and_check_weights(weights, X_train_sc)
  pf_list <- lapply(0:max_imp_pow, function(i) 1 / (w^i))
  pf_names <- paste0("1/imp^", 0:max_imp_pow)

  ref_cvm <- NULL
  ref_nz <- NULL
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
    co_list <- coef(cv_best, s = s_choice)
    feat_nonzero <- Reduce(
      "|",
      lapply(co_list, function(cm) as.numeric(cm[-1, 1] != 0))
    )
    n_features <- sum(feat_nonzero)
    coef_obj <- co_list
  }

  list(
    algo = "LLM-Lasso",
    model = best_name,
    method = lam,
    coef = coef_obj,
    n_features = n_features
  )
}


# ------------------------------------------------------------------------------
# BIC-Based Lambda Selection for SSL
#
# Selects the lambda0 value along the SSL solution path that minimizes BIC,
# and returns the corresponding coefficient vector with the intercept prepended.
#
# Arguments:
#   ssl_object - Output of lsp_ssl_map
#   X          - Uncentered predictor matrix
#   y          - Response vector
#
# Returns:
#   Numeric coefficient vector (intercept, beta_1, ..., beta_p)
# ------------------------------------------------------------------------------
select_lambda0_bic <- function(ssl_object, X, y) {
  n <- length(y)
  rss <- apply(ssl_object$beta, 2, function(b) sum((y - X %*% b)^2))
  df <- apply(ssl_object$beta, 2, function(b) sum(b != 0))
  bic <- n * log(rss / n) + df * log(n)

  best_idx <- which.min(bic)
  c(ssl_object$intercept[, best_idx], ssl_object$beta[, best_idx])
}

# ------------------------------------------------------------------------------
# Simulation Functions
#
# All three functions below read the following globals from weight_quality_sims.R:
#   p, s, n, effect_size, y_sd, Xvar, cov_mat, sparsity, a_sigma, b_sigma,
#   tau, iter, burn_in, eta_range, fixed_s, random_s
# ------------------------------------------------------------------------------

# Generate data and fit weight-free baseline models (Lasso, horseshoe, and
# optionally standard SS / SSL without LLM weights). Returns a list containing
# the generated dataset and fitted baseline objects, to be passed to sim_function.
baseline_data_sim_function <- function(seed, n) {
  set.seed(seed)

  X <- MASS::mvrnorm(n, mu = rep(0, p), Xvar * cov_mat)
  beta <- c(rep(0, p - s), rep(effect_size, s))
  alpha <- effect_size
  y <- X %*% beta + alpha + rnorm(n, 0, sd = y_sd)

  lasso_results <- as.vector(coef(glmnet::glmnet(
    X,
    y,
    alpha = 1,
    lambda = glmnet::cv.glmnet(X, y, alpha = 1)$lambda.min
  )))

  hs_fit <- Mhorseshoe::approx_horseshoe(
    y = y,
    X = cbind(1, X),
    burn = 10000,
    iter = 5000
  )
  hs_coef <- hs_fit$BetaHat
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
    data = list(X = X, y = y, alpha = alpha, beta = beta),
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
# simulation replicate. Used to assess sensitivity to the choice of eta.
eta_sensitivity_function <- function(seed, weights, set_eta) {
  set.seed(seed)

  X <- MASS::mvrnorm(n, mu = rep(0, p), Xvar * cov_mat)
  beta <- c(rep(0, p - s), rep(effect_size, s))
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
      weights = weights,
      penalty = "adaptive",
      variance = "fixed"
    ) |>
      select_lambda0_bic(X = X, y = y)
  }

  map(all_fits, ~ compile_model_metrics(.x, weights, alpha, beta))
}
