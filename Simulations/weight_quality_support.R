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

# ex:
# generate_weights(0.9, true_gamma) %>% l1_weight_agreement(true_gamma, .)
# generate_weights(0.9, true_gamma) %>% l2_weight_agreement(true_gamma, .)
# generate_weights(0.9, true_gamma) %>% pairwise_weight_agreement(true_gamma, .)

# helper functions for simulation --------------------------------------------------------

compile_model_metrics <- function(model_object, weights, alpha, beta) {
  signal <- beta != 0

  if (class(model_object) == "numeric") {
    # i.e. lasso objects
    gamma_predict <- as.vector(model_object != 0)[-1]
    l1 <- sum(abs(as.vector(model_object) - c(alpha, beta))) # l1 loss of beta
  } else {
    # mcmc objects
    gamma_predict <- model_object$gamma %>% colMeans()
    l1 <- sum(abs(model_object$beta %>% colMeans() - c(alpha, beta))) # l1 loss of beta
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

# Fixed Sparsity ----------------------------------------------------------------
source("MCMC Samplers/LSP_regression_fixed_s.R")

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
# this function runs the simulation
sim_function <- function(seed, n, weights) {
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

  # run llm-lasso
  llm_lasso_results <- llm_lasso_simp(
    X_train = X,
    y_train = y,
    weights = weights,
    elastic_net = 1,
    regression = TRUE
  )$coef

  metrics_list <- list("lasso" = lasso_results, "llm-lasso" = llm_lasso_results)

  if (fixed_s == TRUE) {
    # run standard discrete spike and slab
    dssr_samples <- fswr_gibbs_sampler(
      X,
      y,
      c = 0,
      sparsity = sparsity,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau
    )

    # run LSP with confidence = 0.5
    rvwr_samples <- fswr_gibbs_sampler(
      X,
      y,
      weights,
      sparsity = sparsity,
      c = weight_confidence,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau
    )

    # run LSP with confidence = 1.0
    cwr_samples <- fswr_gibbs_sampler(
      X,
      y,
      weights,
      sparsity = sparsity,
      c = 1,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau
    )

    metrics_list <- c(
      metrics_list,
      list(
        "standard ss" = dssr_samples,
        "random weights" = rvwr_samples,
        "constant weights" = cwr_samples
      )
    )
  }

  # can also evaluate LSP with random sparsity
  if (random_s == TRUE) {
    # run spike and slab with random s, c = 0
    srsr_samples <- rswr_gibbs_sampler(
      X,
      y,
      weights,
      c = 0,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau,
      a_s = a_s,
      b_s = b_s,
      init_weights = FALSE
    )

    # run LSP with confidence = 0.5
    rvsr_samples <- rswr_gibbs_sampler(
      X,
      y,
      weights,
      c = weight_confidence,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau,
      a_s = a_s,
      b_s = b_s
    )

    # run LSP with confidence = 1.0
    csr_samples <- rswr_gibbs_sampler(
      X,
      y,
      weights,
      c = 1,
      a_sigma = a_sigma,
      b_sigma = b_sigma,
      tau = tau,
      a_s = a_s,
      b_s = b_s
    )

    metrics_list <- c(
      metrics_list,
      list(
        "standard ss, random s" = srsr_samples,
        "random weights, random s" = rvsr_samples,
        "constant weights, random s" = csr_samples
      )
    )
  }

  # return metrics
  metrics <- map(
    metrics_list,
    ~ compile_model_metrics(.x, weights, alpha, beta)
  )

  return(metrics)
}
