# AKI Analysis — Model Fitting and Evaluation Support
#
# Utility functions for the AKI real-data analysis. Covers data preparation,
# cross-validation partitioning, and six train-and-evaluate routines spanning
# different model families and weight integration strategies:
#   train_and_evaluate_baselines           — weight-free baseline models
#   train_and_evaluate_fixed_eta           — LSP at a user-specified eta
#   train_and_evaluate_random_eta          — LSP with discrete uniform prior on eta
#   train_and_evaluate_probability_weights — LLM weights used as direct inclusion probs
#   train_and_evaluate_non_ss              — SSL and LLM-Lasso comparisons
#   train_and_evaluate_spike_and_slab      — SS comparisons across weight strategies

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library("dplyr")
  library("purrr")
  library("tibble")
  library("readr")
  library("rsample")
  library("furrr")
  library("future")
})

`%!in%` <- Negate(`%in%`)

source("LSP_SS/LSP_SSR_fixed_s.R")
source("LSP_SS/LSP_SSR_random_s.R")
source("LSP_SSL/LSP_SSLR.R")


# Process Weights --------------------------------------------------------

# Retain and order weights to match the columns of a given subgroup dataset
subset_weights <- function(aki_weights_0, aki_data_set, eps = 1e-8) {
  is_prob <- max(aki_weights_0$importance, na.rm = TRUE) <= 1

  out <- aki_weights_0 |>
    select(value, importance) |>
    filter(value %in% colnames(aki_data_set)) |>
    mutate(value = factor(value, levels = colnames(aki_data_set))) |>
    arrange(value)

  if (anyNA(out$importance)) {
    stop(
      "subset_weights: missing importance for ",
      paste(out$value[is.na(out$importance)], collapse = ", ")
    )
  }

  if (is_prob) {
    out <- mutate(out, importance = pmin(pmax(importance, eps), 1 - eps))
  }

  out
}

# ------------------------------------------------------------------------------
# LLM-Lasso
#
# Adapted from Zhang et al.: https://github.com/pilancilab/LLM-Lasso
# Lightly edited for compatibility with this analysis framework.
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

# Area between the candidate CV error curve and the baseline (uniform penalty)
# CV error curve, interpolated to a common sparsity grid. Larger = better.
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

# Fit LLM-Lasso by selecting the penalty factor exponent (1/w^k,
# k = 0,...,max_imp_pow) that maximizes the area between its CV error curve
# and the unweighted baseline, then re-fitting at the chosen penalty.
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

# Fit LLM-Lasso at a single fixed eta value (penalty factor = 1/w^eta) rather
# than searching over a grid of exponents. Useful for sensitivity analysis.
llm_lasso_fixed_eta <- function(
  X_train,
  y_train,
  weights,
  folds_cv = 5,
  elastic_net = 1,
  eta,
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

  cv_best <- glmnet::cv.glmnet(
    x = X_train_sc,
    y = y_train,
    family = glm_family,
    alpha = elastic_net,
    penalty.factor = 1 / (w^eta),
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
    method = lam,
    coef = coef_obj,
    n_features = n_features
  )
}

# ------------------------------------------------------------------------------
# LLM-Select
#
# Adapted from Jeong et al. (2025)
# ------------------------------------------------------------------------------

fit_llm_select <- function(
  X_train,
  y_train,
  weights,
  selection_rate = 0.1
) {
  stopifnot(length(selection_rate) == 1L)

  p <- ncol(X_train)
  k <- max(2L, round(selection_rate * nrow(weights)))

  selected_vars <- weights |>
    slice_max(order_by = importance, n = k, with_ties = FALSE) |>
    pull(value)

  keep <- colnames(X_train) %in% selected_vars
  stopifnot(sum(keep) >= 2L)
  X_sub <- X_train[, keep, drop = FALSE]

  llm_select_fit <- glmnet::glmnet(
    x = X_sub,
    y = y_train,
    alpha = 0,
    lambda = glmnet::cv.glmnet(X_sub, y_train, alpha = 0)$lambda.min
  )

  sub_coef <- as.numeric(coef(llm_select_fit)) # length 1 + sum(keep)
  full_coef <- numeric(p + 1L)
  full_coef[1L] <- sub_coef[1L]
  full_coef[which(keep) + 1L] <- sub_coef[-1L]

  list(
    algo = "LLM-Select",
    selection_rate = selection_rate,
    coef = full_coef
  )
}

# ------------------------------------------------------------------------------
# Data Preparation Utilities
# ------------------------------------------------------------------------------

# Scale a matrix using precomputed training center and scale vectors.
# Missing values (produced by scaling) are imputed with the column mean (0
# on the scaled space).
scale_data <- function(mat, train_center, train_scale) {
  scale(mat, center = train_center, scale = train_scale) %>%
    apply(2, function(x) {
      x[is.na(x)] <- 0
      x
    })
}

# L1 agreement: 1 - mean absolute deviation after min-max scaling
l1_weight_agreement <- function(true_gamma, weights) {
  scaled_weights <- (weights - min(weights)) / (max(weights) - min(weights))
  1 - mean(abs(true_gamma - scaled_weights))
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

# Remove near-constant features: retains only columns where the modal value
# occupies fewer than `threshold` * 100% of observations.
topK_features <- function(data, threshold = 0.9) {
  n <- nrow(data)
  mode_count <- data |> summarize(across(everything(), ~ max(table(.x))))
  data |>
    select(any_of(
      mode_count |> select(where(~ .x < threshold * n)) |> colnames()
    ))
}


# ------------------------------------------------------------------------------
# estimate_phi
#
# Estimates the empirical weight agreement between LLM-derived weights and a
# spike-and-slab posterior inclusion vector fitted on the supplied data. Useful
# as an offline diagnostic for calibrating weight quality before running the
# main analysis.
#
# Arguments:
#   data         - Data frame with outcome and predictors
#   outcome_var  - Name of the outcome column
#   weights      - Data frame with an `importance` column (length p)
#   set_tau      - Slab variance for the spike-and-slab sampler
#   set_sparsity - Prior inclusion probability for the sampler
#   burn_in      - Burn-in iterations
#   iter         - Total sampler iterations
#
# Returns:
#   Named numeric vector with L1 and pairwise weight agreement values
# ------------------------------------------------------------------------------
estimate_phi <- function(
  data,
  outcome_var,
  weights,
  set_tau = 2,
  set_sparsity = 0.05,
  model = "SSL",
  burn_in = 25000,
  iter = 125000
) {
  y_train <- data |> pull(any_of(outcome_var))
  X_train <- data %>% select(-any_of(outcome_var)) %>% as.matrix()

  X_train_scaled <- scale_data(
    X_train,
    colMeans(X_train, na.rm = TRUE),
    apply(X_train, 2, sd, na.rm = TRUE)
  )
  y_train_scaled <- scale_data(
    y_train,
    mean(y_train, na.rm = TRUE),
    sd(y_train, na.rm = TRUE)
  )

  if (model == "SS") {
    ss_results <- lsp_random_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      a_sigma = 1,
      b_sigma = 1,
      tau = set_tau,
      init_weights = FALSE,
      burn_in = burn_in,
      iter = iter,
      E_space = 0
    )

    gamma_est <- as.integer(ss_results$gamma > 0.5)
  } else if (model == "SSL") {
    ssl_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      E_space = 0,
      weights = NULL,
      penalty = "separable",
      variance = "fixed",
      sparsity = set_sparsity
    ) |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)

    gamma_est <- rep(0, ncol(X_train_scaled))
    gamma_est[as.vector(ssl_fit)[-1] > 0] <- 1
  } else if (model == "Lasso") {
    lasso_results <- glmnet::glmnet(
      x = X_train_scaled,
      y = y_train_scaled,
      alpha = 1,
      lambda = glmnet::cv.glmnet(
        X_train_scaled,
        y_train_scaled,
        alpha = 1
      )$lambda.min
    )

    gamma_est <- as.numeric(lasso_results$coef != 0)
  }

  c(
    l1_weight_agreement = l1_weight_agreement(gamma_est, weights$importance),
    ss_pairwise_weight_agreement = pairwise_weight_agreement(
      gamma_est,
      weights$importance
    )
  )
}

# ------------------------------------------------------------------------------
# train_test_split
#
# Partitions data into stratified cross-validation folds and returns a list of
# scaled train/test matrices. When n_folds = 1, all observations are assigned
# to training (assessment set is empty). When n is specified, the training set
# is randomly subsampled to that size before scaling.
#
# Arguments:
#   data        - Data frame with outcome and predictors
#   weights     - Weight data frame to attach to each partition
#   outcome_var - Name of the outcome column
#   n_folds     - Number of CV folds (1 = no held-out test set)
#   repetitions - Number of CV repetitions
#   n           - Optional training set size cap (subsampled if exceeded)
#
# Returns:
#   A list of partition objects, each containing scaled X/y train and test
#   matrices, scale factors for back-transformation, and the weight data frame
# ------------------------------------------------------------------------------
train_test_split <- function(
  data,
  weights,
  outcome_var,
  n_folds = 5,
  repetitions = 1,
  n = NULL
) {
  if (n_folds == 1) {
    split_obj <- make_splits(
      list(analysis = seq_len(nrow(data)), assessment = integer(0)),
      data = data
    )
    folds <- list(splits = list(split_obj))
  } else {
    folds <- vfold_cv(
      data,
      v = n_folds,
      strata = all_of(outcome_var),
      repeats = repetitions
    )
  }

  map(folds$splits, function(x) {
    train_data <- as.data.frame(x)
    test_data <- data[-x[[2]], ]

    if (!is.null(n) && n < nrow(train_data)) {
      train_data <- train_data[sample(nrow(train_data), n), ]
    }

    y_train <- train_data |> pull(any_of(outcome_var))
    X_train <- train_data %>% select(-any_of(outcome_var)) %>% as.matrix()

    X_train_center <- colMeans(X_train, na.rm = TRUE)
    X_train_scale <- apply(X_train, 2, sd, na.rm = TRUE)
    X_train_scaled <- scale_data(X_train, X_train_center, X_train_scale)

    y_train_center <- mean(y_train, na.rm = TRUE)
    y_train_scale <- sd(y_train, na.rm = TRUE)
    y_train_scaled <- scale_data(y_train, y_train_center, y_train_scale)

    X_test_scaled <- test_data |>
      select(-any_of(outcome_var)) |>
      as.matrix() |>
      scale_data(X_train_center, X_train_scale)
    y_test_scaled <- test_data |>
      pull(any_of(outcome_var)) |>
      scale_data(y_train_center, y_train_scale)

    list(
      X_train_scaled = X_train_scaled,
      y_train_scaled = y_train_scaled,
      X_test_scaled = X_test_scaled,
      y_test_scaled = y_test_scaled,
      y_scale_factor = y_train_scale,
      y_loc_factor = y_train_center,
      weights = weights
    )
  })
}

# ------------------------------------------------------------------------------
# Train-and-Evaluate Utilities
# ------------------------------------------------------------------------------

# Per-observation squared error
squared_error <- function(y_true, pred) (y_true - pred)^2

oos_coverage <- function(coverage, level = 0.95) {
  by_fold <- coverage |>
    group_by(method, partition) |>
    summarize(fold_coverage = mean(coverage), .groups = "drop_last") |>
    summarize(mc_se = sd(fold_coverage) / sqrt(n()), .groups = "drop")

  coverage |>
    group_by(method) |>
    summarize(
      n_intervals = n(),
      coverage = mean(coverage),
      mean_width = mean(width),
      median_width = median(width),
      .groups = "drop"
    ) |>
    left_join(by_fold, by = "method") |>
    mutate(nominal = level, .after = method)
}

ss_predictive_summary <- function(
  samples,
  X_test_scaled,
  y_test_scaled,
  y_scale_factor,
  set_tau,
  method,
  level = 0.95
) {
  beta <- samples$beta
  p1 <- ncol(X_test_scaled) + 1L

  is_draws <- is.matrix(beta) && nrow(beta) > 1L && ncol(beta) == p1

  beta_hat <- if (is_draws) as.numeric(colMeans(beta)) else as.numeric(beta)

  if (length(beta_hat) != p1) {
    stop(
      "ss_predictive_summary: beta has length ",
      length(beta_hat),
      " but X_test implies ",
      p1,
      " (method ",
      method,
      ")"
    )
  }

  # can prediction intervals be formed?
  invsigma_2 <- samples$invsigma_2
  n_test <- nrow(X_test_scaled)
  y <- as.numeric(y_test_scaled)

  can_interval <- is_draws &&
    !is.null(invsigma_2) &&
    length(invsigma_2) == nrow(beta) &&
    all(is.finite(invsigma_2)) &&
    all(invsigma_2 > 0) &&
    length(y) == n_test &&
    n_test > 0L

  if (!can_interval) {
    return(list(beta_hat = beta_hat, coverage = tibble()))
  }

  n_draws <- nrow(beta)
  sigma_draws <- sqrt(1 / invsigma_2)

  # rnorm fills column-major, so recycling sigma_draws (length n_draws) gives
  # each draw its own residual sd within every test column.
  pred_mat <- tcrossprod(beta, cbind(1, X_test_scaled)) +
    matrix(
      rnorm(n_draws * n_test, mean = 0, sd = sigma_draws),
      nrow = n_draws,
      ncol = n_test
    )

  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  pred_int <- apply(pred_mat, 2, quantile, probs = probs, names = FALSE)

  list(
    beta_hat = beta_hat,
    coverage = tibble(
      coverage = y >= pred_int[1, ] & y <= pred_int[2, ],
      width = (pred_int[2, ] - pred_int[1, ]) * y_scale_factor,
      method = method
    )
  )
}

# ------------------------------------------------------------------------------
# train_and_evaluate_baselines
#
# Fits weight-free baseline models (Lasso, horseshoe, and optionally standard
# SS / SSL without LLM weights) on a single partition and returns per-
# observation squared errors rescaled to the original response units.
#
# Arguments:
#   seed         - Random seed for reproducibility
#   partition    - Output element from train_test_split
#   fixed_s      - If TRUE, fit SS/SSL with fixed sparsity
#   random_s     - If TRUE, fit SS/SSL with Beta prior on sparsity
#                  (exactly one of fixed_s / random_s should be TRUE)
#   set_tau      - Slab variance
#   set_sparsity - Prior inclusion probability (used when fixed_s = TRUE)
#   set_burn_in  - Burn-in iterations for MCMC samplers
#   set_iter     - Total iterations for MCMC samplers
# ------------------------------------------------------------------------------
train_and_evaluate_baselines <- function(
  seed,
  partition,
  fixed_s = FALSE,
  random_s = TRUE,
  set_tau = 2,
  set_sparsity = 0.05,
  set_burn_in = 25000,
  set_iter = 125000,
  set_thin = 1,
  return_coverage = FALSE,
  set_level = 0.95
) {
  X_train_scaled <- partition$X_train_scaled
  y_train_scaled <- partition$y_train_scaled
  X_test_scaled <- partition$X_test_scaled
  y_test_scaled <- partition$y_test_scaled
  y_scale_factor <- partition$y_scale_factor
  y_loc_factor <- partition$y_loc_factor

  set.seed(seed)

  if (fixed_s) {
    ss_samples <- lsp_fixed_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = NULL,
      sparsity = set_sparsity,
      a_sigma = 1,
      b_sigma = 1,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
      thin = set_thin,
      init_weights = FALSE,
      return_samples = return_coverage
    )
    ssl_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      E_space = 0,
      weights = NULL,
      penalty = "separable",
      variance = "fixed",
      sparsity = set_sparsity
    ) |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)
  } else if (random_s) {
    ss_samples <- lsp_random_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = NULL,
      a_sigma = 1,
      b_sigma = 1,
      a_s = 1,
      b_s = NA,
      s_proposal_sigma = 2,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
      thin = set_thin,
      init_weights = FALSE,
      return_samples = return_coverage
    )
    ssl_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      penalty = "adaptive",
      variance = "fixed",
      E_space = 0,
      weights = NULL
    ) |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)
  }

  ss_summary <- ss_predictive_summary(
    ss_samples,
    X_test_scaled,
    y_test_scaled,
    y_scale_factor,
    set_tau = set_tau,
    method = "ss",
    level = set_level
  )

  lasso_results <- glmnet::glmnet(
    x = X_train_scaled,
    y = y_train_scaled,
    alpha = 1,
    lambda = glmnet::cv.glmnet(
      X_train_scaled,
      y_train_scaled,
      alpha = 1
    )$lambda.min
  )

  hs_fit <- Mhorseshoe::approx_horseshoe(
    y = y_train_scaled,
    X = cbind(1, X_train_scaled),
    burn = 10000,
    iter = 5000
  )
  hs_coef <- hs_fit$BetaHat
  rm(hs_fit)
  gc()

  if (length(y_test_scaled) > 0) {
    mse <- tibble(
      y = (y_test_scaled * y_scale_factor + y_loc_factor)[, 1],
      ss_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% ss_summary$beta_hat
      )[, 1] *
        y_scale_factor^2,
      lasso_squared_error = squared_error(
        y_test_scaled,
        predict(lasso_results, newx = X_test_scaled)
      )[, 1] *
        y_scale_factor^2,
      hs_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% as.matrix(hs_coef)
      )[, 1] *
        y_scale_factor^2,
      ssl_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% ssl_fit
      )[, 1] *
        y_scale_factor^2
    )
  } else {
    mse <- tibble()
  }

  list(mse = mse, coverage = ss_summary$coverage)
}

# ------------------------------------------------------------------------------
# train_and_evaluate_fixed_eta
#
# Fits LSP models (SS and SSL) and LLM-Lasso at a single fixed eta value.
# Used for eta sensitivity analysis.
#
# Arguments:
#   seed         - Random seed
#   partition    - Output element from train_test_split
#   fixed_s      - If TRUE, use fixed-sparsity variants
#   random_s     - If TRUE, use random-sparsity variants
#                  (exactly one of fixed_s / random_s should be TRUE)
#   set_tau      - Slab variance
#   set_eta      - Fixed eta value passed to E_space
#   set_sparsity - Prior inclusion probability (used when fixed_s = TRUE)
#   set_burn_in  - Burn-in iterations
#   set_iter     - Total iterations
# ------------------------------------------------------------------------------
train_and_evaluate_fixed_eta <- function(
  seed,
  partition,
  fixed_s = FALSE,
  random_s = TRUE,
  set_tau = 2,
  set_eta = 1,
  set_sparsity = 0.05,
  set_burn_in = 25000,
  set_iter = 125000,
  set_thin = 1
) {
  X_train_scaled <- partition$X_train_scaled
  y_train_scaled <- partition$y_train_scaled
  X_test_scaled <- partition$X_test_scaled
  y_test_scaled <- partition$y_test_scaled
  weights <- partition$weights
  y_scale_factor <- partition$y_scale_factor
  y_loc_factor <- partition$y_loc_factor

  set.seed(seed)

  llm_lasso_results <- llm_lasso_fixed_eta(
    X_train = X_train_scaled,
    y_train = y_train_scaled,
    weights = weights$importance,
    eta = set_eta,
    elastic_net = 1,
    regression = TRUE
  )$coef

  if (fixed_s) {
    lsp_ss_samples <- lsp_fixed_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = weights$importance,
      E_space = set_eta,
      sparsity = set_sparsity,
      a_sigma = 1,
      b_sigma = 1,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
      thin = set_thin,
      return_samples = FALSE
    )
    lsp_ssl_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      E_space = set_eta,
      weights = weights$importance,
      penalty = "adaptive",
      variance = "fixed",
      sparsity = set_sparsity
    ) |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)
  } else if (random_s) {
    lsp_ss_samples <- lsp_random_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = weights$importance,
      E_space = set_eta,
      a_sigma = 1,
      b_sigma = 1,
      a_s = 1,
      b_s = NA,
      s_proposal_sigma = 2,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
      thin = set_thin,
      return_samples = FALSE
    )
    lsp_ssl_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      E_space = set_eta,
      weights = weights$importance,
      penalty = "adaptive",
      variance = "fixed"
    ) |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)
  }

  if (length(y_test_scaled) > 0) {
    mse <- tibble(
      y = (y_test_scaled * y_scale_factor + y_loc_factor)[, 1],
      lsp_ss_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% lsp_ss_samples$beta
      )[, 1] *
        y_scale_factor^2,
      lsp_ssl_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% lsp_ssl_fit
      )[, 1] *
        y_scale_factor^2,
      llm_lasso_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% llm_lasso_results
      )[, 1] *
        y_scale_factor^2
    )
  } else {
    mse <- tibble()
  }

  list(mse = mse)
}

# ------------------------------------------------------------------------------
# train_and_evaluate_random_eta
#
# Fits LSP models (SS and SSL) and LLM-Lasso with eta selected from a grid.
# This is the primary LSP evaluation function for the main analysis.
#
# Arguments:
#   seed            - Random seed
#   partition       - Output element from train_test_split
#   fixed_s         - If TRUE, use fixed-sparsity variants
#   random_s        - If TRUE, use random-sparsity variants
#                     (exactly one of fixed_s / random_s should be TRUE)
#   set_tau         - Slab variance
#   set_eta_range   - Grid of eta values passed to E_space
#   set_sparsity    - Prior inclusion probability (used when fixed_s = TRUE)
#   set_burn_in     - Burn-in iterations
#   set_iter        - Total iterations
# ------------------------------------------------------------------------------
train_and_evaluate_random_eta <- function(
  seed,
  partition,
  random_s = TRUE,
  fixed_s = FALSE,
  set_tau = 1,
  set_eta_range = 1,
  set_sparsity = 0.05,
  set_burn_in = 25000,
  set_iter = 125000,
  set_thin = 1,
  return_coverage = FALSE,
  set_level = 0.95
) {
  X_train_scaled <- partition$X_train_scaled
  y_train_scaled <- partition$y_train_scaled
  X_test_scaled <- partition$X_test_scaled
  y_test_scaled <- partition$y_test_scaled
  weights <- partition$weights
  y_scale_factor <- partition$y_scale_factor
  y_loc_factor <- partition$y_loc_factor

  set.seed(seed)

  llm_lasso_results <- llm_lasso_simp(
    X_train = X_train_scaled,
    y_train = y_train_scaled,
    weights = weights$importance,
    elastic_net = 1,
    regression = TRUE
  )$coef

  if (fixed_s) {
    lsp_ss_samples <- lsp_fixed_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = weights$importance,
      E_space = set_eta_range,
      sparsity = set_sparsity,
      a_sigma = 1,
      b_sigma = 1,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
      thin = set_thin,
      return_samples = return_coverage
    )
    lsp_ssl_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      E_space = set_eta_range,
      weights = weights$importance,
      penalty = "separable",
      variance = "fixed",
      sparsity = set_sparsity
    ) |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)
  } else if (random_s) {
    lsp_ss_samples <- lsp_random_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = weights$importance,
      E_space = set_eta_range,
      a_sigma = 1,
      b_sigma = 1,
      a_s = 1,
      b_s = NA,
      s_proposal_sigma = 2,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
      thin = set_thin,
      return_samples = return_coverage
    )
    lsp_ssl_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      E_space = set_eta_range,
      weights = weights$importance,
      penalty = "adaptive",
      variance = "fixed"
    ) |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)
  }

  lsp_ss_summary <- ss_predictive_summary(
    lsp_ss_samples,
    X_test_scaled,
    y_test_scaled,
    y_scale_factor,
    set_tau = set_tau,
    method = "lsp_ss",
    level = set_level
  )

  if (length(y_test_scaled) > 0) {
    mse <- tibble(
      y = (y_test_scaled * y_scale_factor + y_loc_factor)[, 1],
      lsp_ss_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% lsp_ss_summary$beta_hat
      )[, 1] *
        y_scale_factor^2,
      lsp_ssl_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% lsp_ssl_fit
      )[, 1] *
        y_scale_factor^2,
      llm_lasso_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% llm_lasso_results
      )[, 1] *
        y_scale_factor^2
    )
  } else {
    mse <- tibble()
  }

  list(mse = mse, coverage = lsp_ss_summary$coverage)
}

# ------------------------------------------------------------------------------
# train_and_evaluate_probability_weights
#
# Fits SS and SSL models in which the LLM importance weights are used directly
# as per-covariate prior inclusion probabilities.
# ------------------------------------------------------------------------------
train_and_evaluate_probability_weights <- function(
  seed,
  partition,
  set_tau = 1,
  set_burn_in = 25000,
  set_iter = 125000,
  set_thin = 1,
  return_coverage = FALSE,
  set_level = 0.95
) {
  X_train_scaled <- partition$X_train_scaled
  y_train_scaled <- partition$y_train_scaled
  X_test_scaled <- partition$X_test_scaled
  y_test_scaled <- partition$y_test_scaled
  weights <- partition$weights
  y_scale_factor <- partition$y_scale_factor
  y_loc_factor <- partition$y_loc_factor

  set.seed(seed)

  # Weights supplied as sparsity: each theta_j = importance_j (no eta scaling)
  ss_prob_samples <- lsp_fixed_ss_gibbs_sampler(
    X_train_scaled,
    y_train_scaled,
    weights = NULL,
    E_space = 0,
    sparsity = weights$importance,
    a_sigma = 1,
    b_sigma = 1,
    tau = set_tau,
    burn_in = set_burn_in,
    iter = set_iter,
    thin = set_thin,
    return_samples = return_coverage
  )

  ssl_prob_fit <- lsp_ssl_map(
    X_train_scaled,
    y_train_scaled,
    penalty = "separable",
    variance = "fixed",
    E_space = 0,
    weights = NULL,
    sparsity = weights$importance
  ) |>
    select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)

  ss_prob_summary <- ss_predictive_summary(
    ss_prob_samples,
    X_test_scaled,
    y_test_scaled,
    y_scale_factor,
    set_tau = set_tau,
    method = "ss_prob",
    level = set_level
  )

  if (length(y_test_scaled) > 0) {
    mse <- tibble(
      y = (y_test_scaled * y_scale_factor + y_loc_factor)[, 1],
      ss_prob_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% ss_prob_summary$beta_hat
      )[, 1] *
        y_scale_factor^2,
      ssl_prob_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% ssl_prob_fit
      )[, 1] *
        y_scale_factor^2
    )
  } else {
    mse <- tibble()
  }

  list(mse = mse, coverage = ss_prob_summary$coverage)
}

# ------------------------------------------------------------------------------
# train_and_evaluate_non_ss
#
# Fits non-MCMC comparators: Lasso, horseshoe, LLM-Lasso, LLM-Select and SSL variants
# (with and without LLM weights, using both BIC and a fixed lambda0 = n/2
# selection rule).
#
# Arguments:
#   seed           - Random seed
#   partition      - Partition from train_test_split (discretized weights)
#   prob_partition - Partition from train_test_split (direct probability weights)
#   fixed_s        - If TRUE, fit fixed-sparsity SSL variants
#   random_s       - If TRUE, fit random-sparsity SSL variants
#   set_eta_range  - Eta grid for LSP-SSL
#   set_sparsity   - Prior inclusion probability (used when fixed_s = TRUE)
# ------------------------------------------------------------------------------
train_and_evaluate_non_ss <- function(
  seed,
  partition,
  prob_partition,
  fixed_s = FALSE,
  random_s = TRUE,
  set_eta_range = NULL,
  set_sparsity = 0.05,
  llm_select_selection_rate = c(0.10, 0.30)
) {
  X_train_scaled <- partition$X_train_scaled
  y_train_scaled <- partition$y_train_scaled
  X_test_scaled <- partition$X_test_scaled
  y_test_scaled <- partition$y_test_scaled
  weights <- partition$weights
  prob_weights <- prob_partition$weights
  y_scale_factor <- partition$y_scale_factor
  y_loc_factor <- partition$y_loc_factor

  set.seed(seed)

  lasso_results <- glmnet::glmnet(
    x = X_train_scaled,
    y = y_train_scaled,
    alpha = 1,
    lambda = glmnet::cv.glmnet(
      X_train_scaled,
      y_train_scaled,
      alpha = 1
    )$lambda.min
  )

  hs_fit <- Mhorseshoe::approx_horseshoe(
    y = y_train_scaled,
    X = cbind(1, X_train_scaled),
    burn = 10000,
    iter = 5000
  )
  hs_coef <- hs_fit$BetaHat
  rm(hs_fit)
  gc()

  llm_lasso_results <- llm_lasso_simp(
    X_train = X_train_scaled,
    y_train = y_train_scaled,
    weights = weights$importance,
    elastic_net = 1,
    regression = TRUE
  )$coef

  llm_select_results <- map(llm_select_selection_rate, function(sr) {
    fit_llm_select(
      X_train = X_train_scaled,
      y_train = y_train_scaled,
      weights = weights,
      selection_rate = sr
    )$coef
  }) |>
    set_names(paste0(
      "llm_select_",
      round(llm_select_selection_rate * 100, digits = 2),
      "_squared_error"
    ))

  if (fixed_s) {
    ssl_fixed_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      E_space = 0,
      weights = NULL,
      penalty = "separable",
      variance = "fixed",
      sparsity = set_sparsity
    )
    ssl_bic_fixed <- ssl_fixed_fit |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)
    ssl_n_2_fixed <- c(ssl_fixed_fit$intercept[50], ssl_fixed_fit$beta[, 50])

    lsp_ssl_fixed_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      E_space = set_eta_range,
      weights = weights$importance,
      penalty = "separable",
      variance = "fixed",
      sparsity = set_sparsity
    )
    lsp_ssl_bic_fixed <- lsp_ssl_fixed_fit |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)
    lsp_ssl_n_2_fixed <- c(
      lsp_ssl_fixed_fit$intercept[50],
      lsp_ssl_fixed_fit$beta[, 50]
    )
  }

  if (random_s) {
    ssl_random_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      penalty = "adaptive",
      variance = "fixed",
      E_space = 0,
      weights = NULL
    )
    ssl_bic_random <- ssl_random_fit |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)
    ssl_n_2_random <- c(ssl_random_fit$intercept[50], ssl_random_fit$beta[, 50])

    lsp_ssl_random_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      E_space = set_eta_range,
      weights = weights$importance,
      penalty = "adaptive",
      variance = "fixed"
    )
    lsp_ssl_bic_random <- lsp_ssl_random_fit |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)
    lsp_ssl_n_2_random <- c(
      lsp_ssl_random_fit$intercept[50],
      lsp_ssl_random_fit$beta[, 50]
    )
  }

  # Probability-weight SSL: importance weights used directly as sparsity
  ssl_prob_fit <- lsp_ssl_map(
    X_train_scaled,
    y_train_scaled,
    penalty = "separable",
    variance = "fixed",
    E_space = 0,
    weights = NULL,
    sparsity = prob_weights$importance
  )
  ssl_prob_bic <- ssl_prob_fit |>
    select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)
  ssl_prob_n_2 <- c(
    ssl_prob_fit$intercept[50],
    ssl_prob_fit$beta[, 50]
  )

  if (length(y_test_scaled) > 0) {
    se <- function(coef_vec) {
      squared_error(y_test_scaled, cbind(1, X_test_scaled) %*% coef_vec)[, 1] *
        y_scale_factor^2
    }

    mse_list <- list(
      y = (y_test_scaled * y_scale_factor + y_loc_factor)[, 1],
      lasso_squared_error = squared_error(
        y_test_scaled,
        predict(lasso_results, newx = X_test_scaled)
      )[, 1] *
        y_scale_factor^2,
      hs_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% as.matrix(hs_coef)
      )[, 1] *
        y_scale_factor^2,
      llm_lasso_squared_error = se(llm_lasso_results),
      ssl_prob_bic_squared_error = se(ssl_prob_bic),
      ssl_prob_n_2_squared_error = se(ssl_prob_n_2)
    ) |>
      c(map(llm_select_results, se))

    if (fixed_s) {
      mse_list$ssl_bic_fixed_squared_error <- se(ssl_bic_fixed)
      mse_list$ssl_n_2_fixed_squared_error <- se(ssl_n_2_fixed)
      mse_list$lsp_ssl_bic_fixed_squared_error <- se(lsp_ssl_bic_fixed)
      mse_list$lsp_ssl_n_2_fixed_squared_error <- se(lsp_ssl_n_2_fixed)
    }
    if (random_s) {
      mse_list$ssl_bic_random_squared_error <- se(ssl_bic_random)
      mse_list$ssl_n_2_random_squared_error <- se(ssl_n_2_random)
      mse_list$lsp_ssl_bic_random_squared_error <- se(lsp_ssl_bic_random)
      mse_list$lsp_ssl_n_2_random_squared_error <- se(lsp_ssl_n_2_random)
    }

    mse <- tibble::as_tibble(mse_list)
  } else {
    mse <- tibble()
  }

  list(mse = mse)
}

# ------------------------------------------------------------------------------
# train_and_evaluate_spike_and_slab
#
# Fits and compares three SS variants on a single partition:
#   (1) Standard SS (no weights)
#   (2) LSP-SS (LLM weights, eta prior)
#   (3) Probability-weight SS (importance weights as direct inclusion probs)
#
# Arguments:
#   seed           - Random seed
#   partition      - Partition from train_test_split (discretized weights)
#   prob_partition - Partition from train_test_split (probability importance weights)
#   fixed_s        - If TRUE, use fixed-sparsity samplers
#   random_s       - If TRUE, use random-sparsity samplers
#                    (exactly one of fixed_s / random_s should be TRUE)
#   set_tau        - Slab variance
#   set_eta_range  - Eta grid for LSP-SS
#   set_sparsity   - Prior inclusion probability (used when fixed_s = TRUE)
#   set_burn_in    - Burn-in iterations
#   set_iter       - Total iterations
# ------------------------------------------------------------------------------
train_and_evaluate_spike_and_slab <- function(
  seed,
  partition,
  prob_partition,
  random_s = TRUE,
  fixed_s = FALSE,
  set_tau = 2,
  set_eta_range = NULL,
  set_sparsity = 0.05,
  set_burn_in = 25000,
  set_iter = 125000,
  set_thin = 1,
  return_coverage = FALSE,
  set_level = 0.95
) {
  X_train_scaled <- partition$X_train_scaled
  y_train_scaled <- partition$y_train_scaled
  X_test_scaled <- partition$X_test_scaled
  y_test_scaled <- partition$y_test_scaled
  weights <- partition$weights
  prob_weights <- prob_partition$weights
  y_scale_factor <- partition$y_scale_factor
  y_loc_factor <- partition$y_loc_factor

  set.seed(seed)

  if (fixed_s) {
    ss_samples <- lsp_fixed_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = NULL,
      sparsity = set_sparsity,
      a_sigma = 1,
      b_sigma = 1,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
      thin = set_thin,
      init_weights = FALSE,
      return_samples = return_coverage
    )
    lsp_ss_samples <- lsp_fixed_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = weights$importance,
      E_space = set_eta_range,
      sparsity = set_sparsity,
      a_sigma = 1,
      b_sigma = 1,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
      thin = set_thin,
      return_samples = return_coverage
    )
  } else if (random_s) {
    ss_samples <- lsp_random_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = NULL,
      a_sigma = 1,
      b_sigma = 1,
      a_s = 1,
      b_s = NA,
      s_proposal_sigma = 2,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
      thin = set_thin,
      init_weights = FALSE,
      return_samples = return_coverage
    )
    lsp_ss_samples <- lsp_random_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = weights$importance,
      E_space = set_eta_range,
      a_sigma = 1,
      b_sigma = 1,
      a_s = 1,
      b_s = NA,
      s_proposal_sigma = 2,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
      thin = set_thin,
      return_samples = return_coverage
    )
  }

  # Probability-weight SS: importance weights used directly as sparsity
  ss_prob_samples <- lsp_fixed_ss_gibbs_sampler(
    X_train_scaled,
    y_train_scaled,
    weights = NULL,
    E_space = 0,
    sparsity = prob_weights$importance,
    a_sigma = 1,
    b_sigma = 1,
    tau = set_tau,
    burn_in = set_burn_in,
    iter = set_iter,
    thin = set_thin,
    return_samples = return_coverage
  )

  ss_summary <- ss_predictive_summary(
    ss_samples,
    X_test_scaled,
    y_test_scaled,
    y_scale_factor,
    set_tau = set_tau,
    method = "ss",
    level = set_level
  )

  lsp_ss_summary <- ss_predictive_summary(
    lsp_ss_samples,
    X_test_scaled,
    y_test_scaled,
    y_scale_factor,
    set_tau = set_tau,
    method = "lsp_ss",
    level = set_level
  )

  ss_prob_summary <- ss_predictive_summary(
    ss_prob_samples,
    X_test_scaled,
    y_test_scaled,
    y_scale_factor,
    set_tau = set_tau,
    method = "ss_prob",
    level = set_level
  )

  if (length(y_test_scaled) > 0) {
    mse <- tibble(
      y = (y_test_scaled * y_scale_factor + y_loc_factor)[, 1],
      ss_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% ss_summary$beta_hat
      )[, 1] *
        y_scale_factor^2,
      lsp_ss_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% lsp_ss_summary$beta_hat
      )[, 1] *
        y_scale_factor^2,
      ss_prob_squared_error = squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% ss_prob_summary$beta_hat
      )[, 1] *
        y_scale_factor^2
    )
  } else {
    mse <- tibble()
  }

  list(
    mse = mse,
    coverage = bind_rows(
      ss_summary$coverage,
      lsp_ss_summary$coverage,
      ss_prob_summary$coverage
    )
  )
}
