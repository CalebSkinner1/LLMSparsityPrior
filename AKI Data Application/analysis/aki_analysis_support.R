# aki_analysis_support

# load libraries
suppressPackageStartupMessages({
  library("tidyverse")
  library("tidymodels")
  library("furrr")
  library("future")
})

# negate function
`%!in%` <- Negate(`%in%`)

# load Spike-and-Slab Functions
source("LSP_SS/LSP_SSR_fixed_s.R")
source("LSP_SS/LSP_SSR_random_s.R")

# load Spike-and-Slab LASSO Functions
source("LSP_SSL/LSP_SSLR.R")

# LLM-Lasso ---------------------------------------------------------------
# code may be found at https://github.com/pilancilab/LLM-Lasso, lightly edited for functionality

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
    co_list <- coef(cv_best, s = s_choice) # list (one coef matrix per class)
    nz_by_class <- lapply(co_list, function(cm) as.numeric(cm[-1, 1] != 0))
    feat_nonzero <- Reduce("|", nz_by_class)
    n_features <- sum(feat_nonzero)
    coef_obj <- co_list
  }

  list(
    algo = "LLM-Lasso",
    method = lam,
    coef = coef_obj, # numeric vector (gaussian/binomial) or list per class (multinomial)
    n_features = n_features
  )
}

# Data Prep ---------------------------------------------------------------

# scaling/imputation function
scale_data <- function(mat, train_center, train_scale) {
  scale(mat, center = train_center, scale = train_scale) %>% # scale data
    apply(., 2, function(x) {
      # replace NAs with mean
      x[is.na(x)] <- 0
      return(x)
    })
}

# compute l1 weight agreement
l1_weight_agreement <- function(true_gamma, weights) {
  # scale weights
  scaled_weights <- (weights - min(weights)) / (max(weights) - min(weights))

  # l1 distance
  l1 <- mean(abs(true_gamma - scaled_weights))

  # rescale entropy on 0 to 1 scale: "agreement"
  1 - l1
}

# compute pairwise weight agreement
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

# Select top 200 correlative features to remove degenerate columns -------
# function: subset data
topK_features <- function(data, threshold = 0.9) {
  # observations
  n <- data |> nrow()

  # number of times the mode occurs in each feature
  mode_count <- data |>
    summarize(
      across(everything(), ~ max(table(.x)))
    )

  # remove features that do not display unique values at least threshold percent of the time
  selected_features <- mode_count |>
    select(where(~ .x < threshold * n)) |>
    colnames()

  # return data
  data |> select(any_of(selected_features))
}

# estimate phi
estimate_phi <- function(
  data,
  outcome_var,
  weights,
  set_tau = 2,
  set_sparsity = 0.05,
  burn_in = 25000,
  iter = 125000
) {
  # create training data
  y_train <- data |> pull(any_of(outcome_var))
  X_train <- data %>% select(-any_of(outcome_var)) %>% as.matrix()

  X_train_center <- X_train %>% colMeans(na.rm = TRUE)
  X_train_scale <- X_train %>% apply(., 2, sd, na.rm = TRUE)
  X_train_scaled <- scale_data(X_train, X_train_center, X_train_scale)

  y_train_center <- y_train |> mean(na.rm = TRUE)
  y_train_scale <- y_train |> sd(na.rm = TRUE)
  y_train_scaled <- scale_data(y_train, y_train_center, y_train_scale)

  # lasso_results <- glmnet::glmnet(
  #   x = X_train_scaled,
  #   y = y_train_scaled,
  #   alpha = 1,
  #   lambda = (glmnet::cv.glmnet(
  #     X_train_scaled,
  #     y_train_scaled,
  #     alpha = 1
  #   )$lambda.min)
  # )

  # lasso_gamma <- rep(0, nrow(weights))
  # lasso_gamma[lasso_results$beta[, 1] != 0] <- 1

  # lasso_weight_agreement_l1 <- l1_weight_agreement(lasso_gamma, weights$importance)

  ss_results <- lsp_fixed_ss_gibbs_sampler(
    X_train_scaled,
    y_train_scaled,
    sparsity = set_sparsity,
    a_sigma = 1,
    b_sigma = 1,
    tau = set_tau,
    init_weight = FALSE,
    burn_in = burn_in,
    iter = iter,
    c = 0
  )

  ss_gamma <- rep(0, nrow(weights))
  ss_gamma[ss_results$gamma > 0.5] <- 1

  ss_weight_agreement_l1 <- l1_weight_agreement(ss_gamma, weights$importance)
  ss_weight_agreement_pairwise <- pairwise_weight_agreement(
    ss_gamma,
    weights$importance
  )

  c(
    "ss_l1_weight_agreement" = ss_weight_agreement_l1,
    "ss_pairwise_weight_agreement" = ss_weight_agreement_pairwise
  )
}

# train-test split
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
      list(
        analysis = seq_len(nrow(data)),
        assessment = integer(0)
      ),
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

  centered_partitions <- map(
    folds$splits,
    function(x) {
      train_data <- as.data.frame(x)
      test_data <- data[-x[[2]], ]

      # subsample training set if n is specified
      if (!is.null(n) && n < nrow(train_data)) {
        train_data <- train_data[sample(nrow(train_data), n), ]
      }

      # create training data
      y_train <- train_data |> pull(any_of(outcome_var))
      X_train <- train_data %>% select(-any_of(outcome_var)) %>% as.matrix()

      X_train_center <- X_train %>% colMeans(na.rm = TRUE)
      X_train_scale <- X_train %>% apply(., 2, sd, na.rm = TRUE)
      X_train_scaled <- scale_data(X_train, X_train_center, X_train_scale)

      y_train_center <- y_train |> mean(na.rm = TRUE)
      y_train_scale <- y_train |> sd(na.rm = TRUE)
      y_train_scaled <- scale_data(y_train, y_train_center, y_train_scale)

      # create testing data
      X_test_scaled <- test_data |>
        select(-any_of(outcome_var)) |>
        as.matrix() |>
        scale_data(X_train_center, X_train_scale)
      y_test_scaled <- test_data |>
        pull(any_of(outcome_var)) |>
        scale_data(y_train_center, y_train_scale)

      list(
        "X_train_scaled" = X_train_scaled,
        "y_train_scaled" = y_train_scaled,
        "X_test_scaled" = X_test_scaled,
        "y_test_scaled" = y_test_scaled,
        "y_scale_factor" = y_train_scale,
        "y_loc_factor" = y_train_center,
        "weights" = weights
      )
    }
  )
  return(centered_partitions)
}

# Train and Evaluate Models ----------------------------------------------
squared_error <- function(y_true, pred) {
  # squared error
  (y_true - pred)^2
}

select_lambda0_bic <- function(ssl_object, X, y) {
  n <- length(y)

  rss <- apply(
    ssl_object$beta,
    2,
    function(b) sum((y - X %*% b)^2)
  )
  df <- apply(
    ssl_object$beta,
    2,
    function(b) sum(b != 0)
  )
  bic <- n * log(rss / n) + df * log(n)
  best_idx <- which.min(bic)

  c(ssl_object$intercept[, best_idx], ssl_object$beta[, best_idx])
}

train_and_evaluate_baselines <- function(
  seed,
  partition,
  fixed_s = FALSE,
  random_s = TRUE,
  set_tau = 2,
  set_sparsity = 0.05,
  set_burn_in = 25000,
  set_iter = 125000
) {
  X_train_scaled <- partition$X_train_scaled
  y_train_scaled <- partition$y_train_scaled
  X_test_scaled <- partition$X_test_scaled
  y_test_scaled <- partition$y_test_scaled
  weights <- partition$weights
  y_scale_factor <- partition$y_scale_factor
  y_loc_factor <- partition$y_loc_factor

  set.seed(seed)
  if (fixed_s) {
    # spike and slab
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
      init_weights = FALSE,
      return_samples = FALSE
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
  }
  if (random_s) {
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
      init_weights = FALSE,
      return_samples = FALSE
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

  # lasso
  lasso_results <- glmnet::glmnet(
    x = X_train_scaled,
    y = y_train_scaled,
    alpha = 1,
    lambda = (glmnet::cv.glmnet(
      X_train_scaled,
      y_train_scaled,
      alpha = 1
    )$lambda.min)
  )

  # horseshoe prior
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
    ss_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% ss_samples$beta
    )
    lasso_se <- squared_error(
      y_test_scaled,
      predict(lasso_results, newx = X_test_scaled)
    )
    hs_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% as.matrix(hs_coef)
    )
    ssl_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% ssl_fit
    )

    # store mse
    mse <- tibble(
      y = (y_test_scaled * y_scale_factor + y_loc_factor)[, 1],
      ss_squared_error = ss_se[, 1] * (y_scale_factor^2),
      lasso_squared_error = lasso_se[, 1] * (y_scale_factor^2),
      hs_squared_error = hs_se[, 1] * (y_scale_factor^2),
      ssl_squared_error = ssl_se[, 1] * (y_scale_factor^2)
    )
  } else {
    mse <- tibble()
  }

  list(
    "mse" = mse
  )
}

train_and_evaluate_fixed_eta <- function(
  seed,
  partition,
  fixed_s = FALSE,
  random_s = TRUE,
  set_tau = 2,
  set_eta = 1,
  set_sparsity = 0.05,
  set_confidence = 1,
  set_burn_in = 25000,
  set_iter = 125000
) {
  X_train_scaled <- partition$X_train_scaled
  y_train_scaled <- partition$y_train_scaled
  X_test_scaled <- partition$X_test_scaled
  y_test_scaled <- partition$y_test_scaled
  weights <- partition$weights
  y_scale_factor <- partition$y_scale_factor
  y_loc_factor <- partition$y_loc_factor

  set.seed(seed)
  # llm-lasso
  llm_lasso_results <- llm_lasso_fixed_eta(
    X_train = X_train_scaled,
    y_train = y_train_scaled,
    weights = weights$importance,
    eta = set_eta,
    elastic_net = 1,
    regression = TRUE
  )$coef

  if (fixed_s) {
    # spike and slab constant weights
    lsp_ss_samples <- lsp_fixed_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = weights$importance,
      c = set_confidence,
      eta = set_eta,
      s = set_sparsity,
      a_sigma = 1,
      b_sigma = 1,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
      return_samples = FALSE
    )

    lsp_ssl_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      E_space = set_eta,
      weights = weights$importance,
      penalty = "separable",
      variance = "fixed",
      sparsity = set_sparsity
    ) |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)
  }
  if (random_s) {
    # spike and slab constant weights
    lsp_ss_samples <- lsp_random_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = weights$importance,
      c = set_confidence,
      eta = set_eta,
      a_sigma = 1,
      b_sigma = 1,
      a_s = 1,
      b_s = NA,
      s_proposal_sigma = 2,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
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
    lsp_ss_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% lsp_ss_samples$beta
    )

    lsp_ssl_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% lsp_ssl_fit
    )

    llm_lasso_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% llm_lasso_results
    )

    # store mse
    mse <- tibble(
      y = (y_test_scaled * y_scale_factor + y_loc_factor)[, 1],
      # rescaling the mse by the original units
      lsp_ss_squared_error = lsp_ss_se[, 1] * (y_scale_factor^2),
      lsp_ssl_squared_error = lsp_ssl_se[, 1] * (y_scale_factor^2),
      llm_lasso_squared_error = llm_lasso_se[, 1] * (y_scale_factor^2),
    )
  } else {
    mse <- tibble()
  }

  list(
    "mse" = mse
  )
}

train_and_evaluate_random_eta <- function(
  seed,
  partition,
  random_s = TRUE,
  fixed_s = FALSE,
  set_tau = 1,
  set_eta_range = 1,
  set_sparsity = 0.05,
  set_confidence_range = 1,
  set_burn_in = 25000,
  set_iter = 125000
) {
  X_train_scaled <- partition$X_train_scaled
  y_train_scaled <- partition$y_train_scaled
  X_test_scaled <- partition$X_test_scaled
  y_test_scaled <- partition$y_test_scaled
  weights <- partition$weights
  y_scale_factor <- partition$y_scale_factor
  y_loc_factor <- partition$y_loc_factor

  set.seed(seed)

  # llm-lasso
  llm_lasso_results <- llm_lasso_simp(
    X_train = X_train_scaled,
    y_train = y_train_scaled,
    weights = weights$importance,
    elastic_net = 1,
    regression = TRUE
  )$coef

  if (fixed_s) {
    # lsp fixed s
    lsp_ss_samples <- lsp_fixed_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = weights$importance,
      E_space = set_eta_range,
      s = set_sparsity,
      a_sigma = 1,
      b_sigma = 1,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
      return_samples = FALSE
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
  }
  if (random_s) {
    # lsp random s
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
      return_samples = FALSE
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
  if (length(y_test_scaled) > 0) {
    lsp_ss_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% lsp_ss_samples$beta
    )
    lsp_ssl_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% lsp_ssl_fit
    )
    llm_lasso_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% llm_lasso_results
    )

    # store mse
    mse <- tibble(
      y = (y_test_scaled * y_scale_factor + y_loc_factor)[, 1],
      # rescaling the mse by the original units
      lsp_ss_squared_error = lsp_ss_se[, 1] * (y_scale_factor^2),
      lsp_ssl_squared_error = lsp_ssl_se[, 1] * (y_scale_factor^2),
      llm_lasso_squared_error = llm_lasso_se[, 1] * (y_scale_factor^2),
    )
  } else {
    mse <- tibble()
  }

  list(
    "mse" = mse
  )
}

train_and_evaluate_continuous_weights <- function(
  seed,
  partition,
  set_tau = 1,
  set_burn_in = 25000,
  set_iter = 125000
) {
  X_train_scaled <- partition$X_train_scaled
  y_train_scaled <- partition$y_train_scaled
  X_test_scaled <- partition$X_test_scaled
  y_test_scaled <- partition$y_test_scaled
  weights <- partition$weights
  y_scale_factor <- partition$y_scale_factor
  y_loc_factor <- partition$y_loc_factor

  set.seed(seed)

  # ss random s
  ss_continuous_samples <- lsp_fixed_ss_gibbs_sampler(
    X_train_scaled,
    y_train_scaled,
    weights = NULL,
    E_space = 0,
    a_sigma = 1,
    b_sigma = 1,
    sparsity = weights$importance, # weights directly set prior inclusion probability
    tau = set_tau,
    burn_in = set_burn_in,
    iter = set_iter,
    return_samples = FALSE
  )

  ssl_continuous_fit <- lsp_ssl_map(
    X_train_scaled,
    y_train_scaled,
    penalty = "separable",
    variance = "fixed",
    E_space = 0,
    weights = NULL,
    sparsity = weights$importance
  ) |>
    select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)

  if (length(y_test_scaled) > 0) {
    ss_cont_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% ss_continuous_samples$beta
    )

    ssl_cont_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% ssl_continuous_fit
    )

    # store mse
    mse <- tibble(
      y = (y_test_scaled * y_scale_factor + y_loc_factor)[, 1],
      # rescaling the mse by the original units
      ss_cont_squared_error = ss_cont_se[, 1] * (y_scale_factor^2),
      ssl_cont_squared_error = ssl_cont_se[, 1] * (y_scale_factor^2)
    )
  } else {
    mse <- tibble()
  }

  list(
    "mse" = mse
  )
}

train_and_evaluate_non_ss <- function(
  seed,
  partition,
  cont_partition,
  fixed_s = FALSE,
  random_s = TRUE,
  set_eta_range = NULL,
  set_sparsity = 0.05
) {
  X_train_scaled <- partition$X_train_scaled
  y_train_scaled <- partition$y_train_scaled
  X_test_scaled <- partition$X_test_scaled
  y_test_scaled <- partition$y_test_scaled
  weights <- partition$weights
  cont_weights <- cont_partition$weights # continuous importance weights
  y_scale_factor <- partition$y_scale_factor
  y_loc_factor <- partition$y_loc_factor

  set.seed(seed)

  # # lasso
  # lasso_results <- glmnet::glmnet(
  #   x = X_train_scaled,
  #   y = y_train_scaled,
  #   alpha = 1,
  #   lambda = (glmnet::cv.glmnet(
  #     X_train_scaled,
  #     y_train_scaled,
  #     alpha = 1
  #   )$lambda.min)
  # )

  # # horseshoe prior
  # hs_fit <- Mhorseshoe::approx_horseshoe(
  #   y = y_train_scaled,
  #   X = cbind(1, X_train_scaled),
  #   burn = 10000,
  #   iter = 5000
  # )
  # hs_coef <- hs_fit$BetaHat
  # rm(hs_fit)
  # gc()

  # # llm-lasso
  # llm_lasso_results <- llm_lasso_simp(
  #   X_train = X_train_scaled,
  #   y_train = y_train_scaled,
  #   weights = weights$importance,
  #   elastic_net = 1,
  #   regression = TRUE
  # )$coef

  if (fixed_s) {
    # spike-and-slab LASSO
    ssl_fixed_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      E_space = 0,
      weights = NULL,
      penalty = "separable",
      variance = "fixed",
      sparsity = set_sparsity
    )

    # select model with bic
    ssl_bic_fixed <- ssl_fixed_fit |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)

    # select model with lambda_0 = n/2 (approximately 50th index)
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
    # select model with bic
    lsp_ssl_bic_fixed <- lsp_ssl_fixed_fit |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)

    # select model with lambda_0 = n/2 (approximately 50th index)
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

    # select model with bic
    ssl_bic_random <- ssl_random_fit |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)

    # select model with lambda_0 = n/2 (approximately 50th index)
    ssl_n_2_random <- c(
      ssl_random_fit$intercept[50],
      ssl_random_fit$beta[, 50]
    )

    lsp_ssl_random_fit <- lsp_ssl_map(
      X_train_scaled,
      y_train_scaled,
      E_space = set_eta_range,
      weights = weights$importance,
      penalty = "adaptive",
      variance = "fixed"
    )
    # select model with bic
    lsp_ssl_bic_random <- lsp_ssl_random_fit |>
      select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)

    # select model with lambda_0 = n/2 (approximately 50th index)
    lsp_ssl_n_2_random <- c(
      lsp_ssl_random_fit$intercept[50],
      lsp_ssl_random_fit$beta[, 50]
    )
  }

  ssl_continuous_fit <- lsp_ssl_map(
    X_train_scaled,
    y_train_scaled,
    penalty = "separable",
    variance = "fixed",
    E_space = 0,
    weights = NULL,
    sparsity = cont_weights$importance
  )

  # select model with bic
  ssl_cont_bic <- ssl_continuous_fit |>
    select_lambda0_bic(X = X_train_scaled, y = y_train_scaled)

  # select model with lambda_0 = n/2 (approximately 50th index)
  ssl_cont_n_2 <- c(
    ssl_continuous_fit$intercept[50],
    ssl_continuous_fit$beta[, 50]
  )

  if (length(y_test_scaled) > 0) {
    # lasso_se <- squared_error(
    #   y_test_scaled,
    #   predict(lasso_results, newx = X_test_scaled)
    # )
    # hs_se <- squared_error(
    #   y_test_scaled,
    #   cbind(1, X_test_scaled) %*% as.matrix(hs_coef)
    # )
    # llm_lasso_se <- squared_error(
    #   y_test_scaled,
    #   cbind(1, X_test_scaled) %*% llm_lasso_results
    # )

    if (fixed_s) {
      ssl_bic_fixed_se <- squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% ssl_bic_fixed
      )
      ssl_n_2_fixed_se <- squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% ssl_n_2_fixed
      )
      lsp_ssl_bic_fixed_se <- squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% lsp_ssl_bic_fixed
      )
      lsp_ssl_n_2_fixed_se <- squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% lsp_ssl_n_2_fixed
      )
    }

    if (random_s) {
      ssl_bic_random_se <- squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% ssl_bic_random
      )
      ssl_n_2_random_se <- squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% ssl_n_2_random
      )
      lsp_ssl_bic_random_se <- squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% lsp_ssl_bic_random
      )
      lsp_ssl_n_2_random_se <- squared_error(
        y_test_scaled,
        cbind(1, X_test_scaled) %*% lsp_ssl_n_2_random
      )
    }
    ssl_cont_bic_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% ssl_cont_bic
    )
    ssl_cont_n_2_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% ssl_cont_n_2
    )

    mse_list <- list(
      y = (y_test_scaled * y_scale_factor + y_loc_factor)[, 1],
      # lasso_squared_error = lasso_se[, 1] * (y_scale_factor^2),
      # hs_squared_error = hs_se[, 1] * (y_scale_factor^2),
      # llm_lasso_squared_error = llm_lasso_se[, 1] * (y_scale_factor^2),
      ssl_cont_bic_squared_error = ssl_cont_bic_se[, 1] * (y_scale_factor^2),
      ssl_cont_n_2_squared_error = ssl_cont_n_2_se[, 1] * (y_scale_factor^2)
    )

    if (fixed_s) {
      mse_list$ssl_bic_fixed_squared_error <- ssl_bic_fixed_se[, 1] *
        (y_scale_factor^2)
      mse_list$ssl_n_2_fixed_squared_error <- ssl_n_2_fixed_se[, 1] *
        (y_scale_factor^2)
      mse_list$lsp_ssl_bic_fixed_squared_error <- lsp_ssl_bic_fixed_se[, 1] *
        (y_scale_factor^2)
      mse_list$lsp_ssl_n_2_fixed_squared_error <- lsp_ssl_n_2_fixed_se[, 1] *
        (y_scale_factor^2)
    }

    if (random_s) {
      mse_list$ssl_bic_random_squared_error <- ssl_bic_random_se[, 1] *
        (y_scale_factor^2)
      mse_list$ssl_n_2_random_squared_error <- ssl_n_2_random_se[, 1] *
        (y_scale_factor^2)
      mse_list$lsp_ssl_bic_random_squared_error <- lsp_ssl_bic_random_se[, 1] *
        (y_scale_factor^2)
      mse_list$lsp_ssl_n_2_random_squared_error <- lsp_ssl_n_2_random_se[, 1] *
        (y_scale_factor^2)
    }

    mse <- tibble::as_tibble(mse_list)
  } else {
    mse <- tibble()
  }

  list("mse" = mse)
}

train_and_evaluate_spike_and_slab <- function(
  seed,
  partition,
  cont_partition,
  random_s = TRUE,
  fixed_s = FALSE,
  set_tau = 2,
  set_eta_range = NULL,
  set_sparsity = 0.05,
  set_burn_in = 25000,
  set_iter = 125000
) {
  X_train_scaled <- partition$X_train_scaled
  y_train_scaled <- partition$y_train_scaled
  X_test_scaled <- partition$X_test_scaled
  y_test_scaled <- partition$y_test_scaled
  weights <- partition$weights
  y_scale_factor <- partition$y_scale_factor
  y_loc_factor <- partition$y_loc_factor
  cont_weights <- cont_partition$weights # continuous importance weights

  set.seed(seed)

  if (fixed_s) {
    # --- Standard Spike-and-Slab (no LLM weights) ---
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
      init_weights = FALSE,
      return_samples = FALSE
    )

    # --- LSP Spike-and-Slab (LLM weights, fixed sparsity) ---
    lsp_ss_samples <- lsp_fixed_ss_gibbs_sampler(
      X_train_scaled,
      y_train_scaled,
      weights = weights$importance,
      E_space = set_eta_range,
      s = set_sparsity,
      a_sigma = 1,
      b_sigma = 1,
      tau = set_tau,
      burn_in = set_burn_in,
      iter = set_iter,
      return_samples = FALSE
    )
  }
  if (random_s) {
    # --- Standard Spike-and-Slab (random sparsity) ---
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
      init_weights = FALSE,
      return_samples = FALSE
    )

    # --- LSP Spike-and-Slab (random sparsity) ---
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
      return_samples = FALSE
    )
  }

  # --- Continuous-weight Spike-and-Slab (always fixed-sparsity style) ---
  ss_continuous_samples <- lsp_fixed_ss_gibbs_sampler(
    X_train_scaled,
    y_train_scaled,
    weights = NULL,
    E_space = 0,
    a_sigma = 1,
    b_sigma = 1,
    sparsity = cont_weights$importance, # weights set inclusion probability directly
    tau = set_tau,
    burn_in = set_burn_in,
    iter = set_iter,
    return_samples = FALSE
  )

  if (length(y_test_scaled) > 0) {
    ss_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% ss_samples$beta
    )
    lsp_ss_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% lsp_ss_samples$beta
    )
    ss_cont_se <- squared_error(
      y_test_scaled,
      cbind(1, X_test_scaled) %*% ss_continuous_samples$beta
    )

    mse <- tibble(
      y = (y_test_scaled * y_scale_factor + y_loc_factor)[, 1],
      ss_squared_error = ss_se[, 1] * (y_scale_factor^2),
      lsp_ss_squared_error = lsp_ss_se[, 1] * (y_scale_factor^2),
      ss_cont_squared_error = ss_cont_se[, 1] * (y_scale_factor^2)
    )
  } else {
    mse <- tibble()
  }

  list("mse" = mse)
}
