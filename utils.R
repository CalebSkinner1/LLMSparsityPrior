# Shared Helper Functions
#
# Helpers used by more than one file in this repository. Sourced by the LSP
# implementations (LSP_SS/, LSP_SSL/), which are in turn sourced by the
# simulation and AKI analysis support files. Provides:
#   - pick_one            -- safe single draw from a vector (LSP SS samplers)
#   - standardize_X       -- column standardization used by all LSP functions
#   - unstandardize_beta  -- maps standardized coefficients to the original scale
#   - LLM-Lasso           -- comparison method (adapted from Zhang et al.) and
#                            its helpers: .scale_like_train,
#                            .align_and_check_weights, cve, llm_lasso_simp

# ------------------------------------------------------------------------------
# Draw one Element of x
#
# Safe replacement for sample(x)[1], which draws from 1:x when x is a single number >= 1.
# ------------------------------------------------------------------------------
pick_one <- function(x) x[sample.int(length(x), 1)]

# ------------------------------------------------------------------------------
# Standardize a Design Matrix
#
# Centers each column of X and scales it to unit (population) variance, so
# that every column has sum of squares n, following the SSLASSO convention.
# Columns with zero variance are centered but left unscaled.
#
# Arguments:
#   X - n x p numeric matrix
#
# Returns:
#   A list with components:
#     X      - Standardized n x p matrix
#     center - Column means of X (length p)
#     scale  - Column scale factors (length p)
# ------------------------------------------------------------------------------
standardize_X <- function(X) {
  center <- colMeans(X)
  X_centered <- sweep(X, 2, center)
  scale <- sqrt(colMeans(X_centered^2))
  scale[!is.finite(scale) | scale == 0] <- 1

  list(X = sweep(X_centered, 2, scale, "/"), center = center, scale = scale)
}

# ------------------------------------------------------------------------------
# Map Standardized Coefficients to the Original Scale
#
# Inverts standardize_X for a coefficient vector (intercept, beta_1, ...,
# beta_p) or for a matrix whose rows are such vectors (e.g., MCMC draws).
#
# Arguments:
#   beta_std - Coefficient vector of length p + 1, or matrix with p + 1 columns,
#              estimated on the standardized design (intercept first)
#   center   - Column means returned by standardize_X
#   scale    - Column scale factors returned by standardize_X
#
# Returns:
#   Coefficients on the original scale of X, in the same shape as beta_std
# ------------------------------------------------------------------------------
unstandardize_beta <- function(beta_std, center, scale) {
  is_vector <- is.null(dim(beta_std))
  beta_mat <- if (is_vector) matrix(beta_std, nrow = 1) else beta_std

  slopes <- sweep(beta_mat[, -1, drop = FALSE], 2, scale, "/")
  intercept <- beta_mat[, 1] - as.vector(slopes %*% center)
  out <- cbind(intercept, slopes, deparse.level = 0)

  if (is_vector) as.vector(out) else out
}

# ------------------------------------------------------------------------------
# LLM-Lasso
#
# Adapted from Zhang et al.: https://github.com/pilancilab/LLM-Lasso
# Lightly edited for compatibility with this repository.
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
  df1 <- tibble(ref_non_zero, ref_cvm) |>
    group_by(ref_non_zero) |>
    summarise(ref_cvm = min(ref_cvm), .groups = "drop") |>
    arrange(ref_non_zero)
  df2 <- tibble(non_zero, cvm) |>
    group_by(non_zero) |>
    summarise(cvm = min(cvm), .groups = "drop") |>
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

  sc <- .scale_like_train(X_train)
  X_train_sc <- sc$X_train
  w <- .align_and_check_weights(weights, X_train_sc)

  # Use one fold assignment for every penalty factor so that the CV error
  # curves are directly comparable
  foldid <- sample(rep_len(seq_len(folds_cv), nrow(X_train_sc)))
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
      foldid = foldid,
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
    foldid = foldid,
    lambda.min.ratio = lambda_min_ratio,
    standardize = FALSE,
    type.measure = type_measure
  )
  s_choice <- if (use_lambda_1se) "lambda.1se" else "lambda.min"
  lam <- if (use_lambda_1se) cv_best$lambda.1se else cv_best$lambda.min

  if (glm_family != "multinomial") {
    co <- as.numeric(coef(cv_best, s = s_choice))
    b <- co[-1] / sc$scale
    b0 <- co[1] - sum(sc$center * b)
    coef_obj <- c(b0, b)
    n_features <- sum(b != 0)
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
