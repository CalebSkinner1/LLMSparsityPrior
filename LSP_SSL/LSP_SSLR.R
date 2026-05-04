# LLM Sparsity Prior (LSP) — Spike-and-Slab Lasso, MAP Estimation
#
# Implements MAP estimation for the LSP-augmented Spike-and-Slab Lasso (SSL)
# via coordinate descent. Adapted from the SSLASSO package (Rockova & George,
# JASA 2018): https://cran.r-project.org/web/packages/SSLASSO/index.html
# Source reference: https://github.com/cran/SSLASSO/blob/master/R/SSLASSO.R
#
# Extensions over the base SSL:
#   - Per-covariate prior inclusion probabilities theta_j = s * v_j, where
#     v_j = w_j^eta / mean(w^eta) are normalized LLM-derived weights.
#   - A Beta(a_s, b_s) prior on global sparsity s, updated within the C routine.
#   - A discrete zero-inflated prior on the concentration parameter eta,
#     selected post-hoc by evaluating the log-posterior across the eta grid.

# ------------------------------------------------------------------------------
# Load compiled C routine
# ------------------------------------------------------------------------------
dll_path <- file.path(getwd(), paste0("LSP_SSL/lsp_ssl", .Platform$dynlib.ext))

if ("lsp_ssl" %in% names(getLoadedDLLs())) {
  dyn.unload(dll_path)
}
dyn.load(dll_path)

stopifnot(is.loaded("SSL_gaussian", PACKAGE = "lsp_ssl"))

# ------------------------------------------------------------------------------
# LSP Spike-and-Slab Lasso — MAP Estimation
#
# Fits the LSP-SSL model along a sequence of lambda0 values (the spike
# penalty). For each lambda0, coordinate descent is run separately for
# each eta in E_space; the eta maximizing the log-posterior is then
# selected, yielding a final coefficient path.
#
# Arguments:
#   X             - n x p predictor matrix (or data frame coercible to matrix)
#   y             - Response vector (length n)
#   weights       - Optional LLM-derived weight vector (length p);
#                   defaults to uniform weights
#   E_space       - Grid of eta values for weight concentration;
#                   NULL triggers automatic grid search
#   eta_zero_mass - Prior mass placed on eta = 0 in the discrete prior on eta;
#                   the remaining mass is split uniformly over eta > 0
#   penalty       - "adaptive" or "separable"; passed to the C routine
#   variance      - "fixed" or "unknown"; governs sigma^2 estimation
#   lambda1       - Slab penalty (scalar); defaults to lambda0[1] if missing
#   lambda0       - Spike penalty sequence (increasing); defaults to a
#                   grid of nlambda values in [1, n]
#   beta.init     - Initial coefficient vector (length p); defaults to zeros
#   nlambda       - Number of lambda0 values when lambda0 is auto-generated
#   sparsity      - Scalar or length-p vector of prior inclusion probabilities;
#                   used as the baseline for theta_j = sparsity_j * v_j
#   sigma         - Initial noise standard deviation (used when variance = "fixed"
#                   or as a starting value when variance = "unknown")
#   a_s           - Shape parameter of the Beta(a_s, b_s) prior on s
#   b_s           - Rate parameter of the Beta(a_s, b_s) prior on s;
#                   defaults to p following Rockova & George (2018)
#   eps           - Convergence tolerance for coordinate descent
#   max.iter      - Maximum coordinate descent iterations per lambda0
#   counter       - Number of EM steps per coordinate descent sweep
#   warn          - If TRUE, emit warnings for non-converged lambda0 values
#
# Returns:
#   An object of class "SSLASSO" (list) with components:
#     beta      - p x nlambda matrix of MAP coefficient estimates
#     intercept - Intercept values along the lambda0 path
#     iter      - Iteration counts at each lambda0
#     lambda0   - Spike penalty sequence used
#     lambda1   - Slab penalty used
#     penalty   - Penalty type used
#     thetas    - Estimated inclusion probabilities along the path
#     sigmas    - Estimated sigma values along the path
#     weights   - Weight vector used
#     E_space   - Eta grid used
#     best_eta  - Selected eta value at each lambda0
#     select    - Binary selection matrix (p x nlambda)
#     model     - Indices of selected variables at the final lambda0
#     n         - Number of observations
# ------------------------------------------------------------------------------
lsp_ssl_map <- function(
  X,
  y,
  weights = NULL,
  E_space = NULL,
  eta_zero_mass = 0.5,
  penalty = c("adaptive", "separable"),
  variance = c("fixed", "unknown"),
  lambda1,
  lambda0,
  beta.init = numeric(ncol(X)),
  nlambda = 100,
  sparsity = 0.01,
  sigma = 1,
  a_s = 1,
  b_s,
  eps = 0.001,
  max.iter = 500,
  counter = 10,
  warn = FALSE
) {
  penalty <- match.arg(penalty)
  variance <- match.arg(variance)

  # --- Input validation and coercion ---

  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix(~ 0 + ., data = X), silent = TRUE)
    if (inherits(tmp, "try-error")) {
      stop("X must be a matrix or able to be coerced to a matrix")
    }
  }
  if (storage.mode(X) == "integer") {
    storage.mode(X) <- "double"
  }

  if (!is.numeric(y)) {
    tmp <- try(y <- as.numeric(y), silent = TRUE)
    if (inherits(tmp, "try-error")) {
      stop("y must be numeric or able to be coerced to numeric")
    }
  }
  if (any(is.na(y)) || any(is.na(X))) {
    stop("Missing data (NA's) detected. Eliminate missing data before calling.")
  }

  XX <- scale(X, center = TRUE, scale = FALSE)
  yy <- y - mean(y)
  p <- ncol(XX)
  n <- length(yy)

  # Validate and expand sparsity to length p
  sparsity <- as.numeric(sparsity)
  if (!length(sparsity) %in% c(1L, p)) {
    stop("sparsity must be length 1 or length ncol(X) (", p, ")")
  }
  if (any(sparsity <= 0) || any(sparsity >= 1)) {
    stop("all sparsity values must be strictly between 0 and 1")
  }
  if (length(sparsity) == 1L) {
    sparsity <- rep(sparsity, p)
  }

  # Validate weights; default to uniform (eta = 0 baseline)
  if (is.null(weights)) {
    weights <- rep(1.0, p)
  } else {
    if (length(weights) != p) {
      stop("weights must have length equal to ncol(X) (", p, ")")
    }
    if (any(weights <= 0)) {
      stop("weights must be strictly positive")
    }
    weights <- as.numeric(weights)
  }

  # --- Build the eta grid ---
  # Search for the largest eta such that all theta_j = sparsity_j * v_j < 1

  if (is.null(E_space)) {
    eta_max <- 0
    step_size <- 1
    theta_bound <- max(sparsity) * max(weights)^eta_max / mean(weights^eta_max)

    while (eta_max <= 20) {
      eta_max <- eta_max + step_size
      theta_bound <- max(sparsity) *
        max(weights)^eta_max /
        mean(weights^eta_max)

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
  } else {
    E_space <- as.numeric(E_space)
    if (length(E_space) < 1) stop("E_space must contain at least one eta value")
  }

  n_eta <- length(E_space)

  # --- Build the lambda0 sequence ---

  if (missing(lambda0)) {
    lambda0 <- seq(1, n, length = nlambda)
    lambda1 <- lambda0[1]
  } else {
    nlambda <- length(lambda0)
    if (missing(lambda1)) lambda1 <- lambda0[1]
  }

  if (sum((lambda0[-1] - lambda0[-nlambda]) > 0) != nlambda - 1) {
    stop("lambda0 must be a monotone increasing sequence")
  }
  if (lambda1 > min(lambda0)) {
    stop("lambda1 must be smaller than lambda0")
  }

  if (missing(b_s)) {
    b_s <- p
  } # Beta(1, p) default: Rockova & George (2018)

  # --- Estimate sigma when variance = "unknown" ---
  # Uses a scaled chi-squared approximation at the sigquant quantile

  df <- 3
  sigquant <- 0.9
  sigest <- sd(yy)
  qchi <- qchisq(1 - sigquant, df)
  ncp <- sigest^2 * qchi / df
  min_sigma2 <- sigest^2 / n
  if (variance == "unknown" && missing(sigma)) {
    sigma <- sqrt(df * ncp / (df + 2))
  }

  # Compute normalized weight vector v for a given eta: v_j = w_j^eta / mean(w^eta)
  make_v_vec <- function(eta) {
    w_eta <- weights^eta
    w_eta / mean(w_eta)
  }

  # --- Run coordinate descent for each eta ---

  eta_results <- vector("list", n_eta)
  for (e in seq_len(n_eta)) {
    v_vec_e <- make_v_vec(E_space[e])
    res_e <- .Call(
      "SSL_gaussian",
      XX,
      yy,
      beta.init,
      penalty,
      variance,
      as.double(lambda1),
      as.numeric(lambda0),
      as.numeric(sparsity),
      as.numeric(v_vec_e),
      as.double(sigma),
      as.double(min_sigma2),
      as.double(a_s),
      as.double(b_s),
      eps,
      as.integer(max.iter),
      as.integer(counter),
      PACKAGE = "lsp_ssl"
    )
    eta_results[[e]] <- list(
      bb = matrix(res_e[[1]], p, nlambda),
      iter = res_e[[3]],
      thetas = res_e[[4]],
      sigmas = res_e[[5]],
      v_vec = v_vec_e
    )
  }

  # --- Select the best eta at each lambda0 via log-posterior comparison ---
  # sigma is held fixed at the eta = 0 estimate to make comparisons across
  # eta values commensurate (only the prior term differs across eta).

  bb <- matrix(0.0, p, nlambda)
  iter <- integer(nlambda)
  thetas <- numeric(nlambda)
  sigmas <- numeric(nlambda)
  best_eta <- numeric(nlambda)

  for (l in seq_len(nlambda)) {
    zero_idx <- which(E_space == 0)
    sigma_ref <- if (length(zero_idx) == 1L) {
      eta_results[[zero_idx]]$sigmas[l]
    } else {
      eta_results[[1]]$sigmas[l]
    }

    # Fall back to the first valid sigma if the reference is missing
    if (is.nan(sigma_ref) || is.na(sigma_ref)) {
      valid <- Filter(
        function(x) !is.nan(x) && !is.na(x),
        lapply(eta_results, function(r) r$sigmas[l])
      )
      sigma_ref <- if (length(valid) > 0) valid[[1]] else 1.0
    }

    log_posts <- vapply(
      seq_len(n_eta),
      function(e) {
        lp <- compute_log_posterior(
          y = yy,
          X = XX,
          beta_final = eta_results[[e]]$bb[, l],
          sigma_final = sigma_ref,
          v_vec = eta_results[[e]]$v_vec,
          lambda1 = lambda1,
          lambda0_final = lambda0[l],
          a_s = a_s,
          b_s = b_s,
          s_fixed = if (penalty == "separable") sparsity else NULL
        )

        # Apply zero-inflated discrete prior on eta: eta = 0 receives mass
        # eta_zero_mass; the remaining mass is split uniformly over eta > 0
        if (n_eta > 1L && E_space[e] != 0.0) {
          lp <- lp +
            log(1 - eta_zero_mass) -
            log(n_eta - 1L) -
            log(eta_zero_mass)
        }
        lp
      },
      numeric(1)
    )

    # If all log-posteriors are NaN (e.g., non-convergence), default to eta = 0
    log_posts[is.nan(log_posts)] <- -Inf
    best_e <- which.max(log_posts)
    if (length(best_e) == 0L) {
      best_e <- 1L
    }

    bb[, l] <- eta_results[[best_e]]$bb[, l]
    iter[l] <- eta_results[[best_e]]$iter[l]
    thetas[l] <- eta_results[[best_e]]$thetas[l]
    sigmas[l] <- eta_results[[best_e]]$sigmas[l]
    best_eta[l] <- E_space[best_e]
  }

  # --- Convergence warnings ---

  if (warn && any(iter == max.iter)) {
    warning(
      "Algorithm did not converge for lambda0: ",
      paste(lambda0[iter == max.iter], collapse = ", ")
    )
  }
  if (iter[nlambda] == max.iter) {
    warning("Algorithm did not converge at the last lambda0 value.")
  }

  # --- Recover intercept and format output ---

  beta <- bb
  intercept <- rep(mean(y), nlambda) - crossprod(attr(XX, "scaled:center"), bb)

  varnames <- if (is.null(colnames(X))) {
    paste0("V", seq_len(ncol(X)))
  } else {
    colnames(X)
  }
  dimnames(beta) <- list(varnames, round(lambda0, digits = 4))

  select <- apply(beta, 2, function(x) as.numeric(x != 0))
  model <- (1:p)[select[, nlambda] == 1]

  structure(
    list(
      beta = beta,
      intercept = intercept,
      iter = iter,
      lambda0 = lambda0,
      lambda1 = lambda1,
      penalty = penalty,
      thetas = thetas,
      sigmas = sigmas,
      weights = weights,
      E_space = E_space,
      best_eta = best_eta,
      select = select,
      model = model,
      n = n
    ),
    class = "SSLASSO"
  )
}


# ------------------------------------------------------------------------------
# Log-Posterior for Spike-and-Slab Lasso
#
# Evaluates the unnormalized log-posterior of the SSL model at the given MAP
# estimates. Used to compare solutions across values of eta.
#
# Arguments:
#   y             - Mean-centered response vector (length n)
#   X             - Mean-centered design matrix (n x p)
#   beta_final    - MAP coefficient vector (length p)
#   sigma_final   - Noise standard deviation at which to evaluate the likelihood
#   v_vec         - Normalized weight vector (length p); v_j = w_j^eta / mean(w^eta)
#   lambda1       - Slab penalty (scalar)
#   lambda0_final - Spike penalty at the current path point (scalar)
#   a_s           - Shape parameter of the Beta prior on s
#   b_s           - Rate parameter of the Beta prior on s
#   s_fixed       - If non-NULL, sparsity is treated as fixed at this value
#                   and the Beta prior on s is excluded (used for "separable" penalty)
#   rss           - Precomputed residual sum of squares; recomputed if NULL
#
# Returns:
#   Scalar log-posterior value
# ------------------------------------------------------------------------------
compute_log_posterior <- function(
  y,
  X,
  beta_final,
  sigma_final,
  v_vec,
  lambda1,
  lambda0_final,
  a_s,
  b_s,
  s_fixed = NULL,
  rss = NULL
) {
  n <- length(y)
  p <- length(beta_final)
  sig2 <- sigma_final^2

  # Sparsity: fixed externally (separable penalty) or estimated from beta
  if (!is.null(s_fixed)) {
    s_val <- s_fixed
    log_prior_s <- 0.0
  } else {
    s_val <- max(
      1e-10,
      min(1.0 - 1e-10, (sum(beta_final != 0.0) + a_s) / (a_s + b_s + p))
    )
    log_prior_s <- (a_s - 1.0) * log(s_val) + (b_s - 1.0) * log(1.0 - s_val)
  }

  theta_v <- pmin(pmax(s_val * v_vec, 1e-8), 1.0 - 1e-8)

  if (is.null(rss)) {
    rss <- sum((y - as.vector(X %*% beta_final))^2)
  }
  log_lik <- -(n / 2.0) * log(sig2) - rss / (2.0 * sig2)

  # Marginal log-prior on beta: mixture of Laplace densities per coefficient
  ab <- abs(beta_final)
  log_prior_beta <- sum(log(pmax(
    theta_v *
      (lambda1 / 2.0) *
      exp(-lambda1 * ab) +
      (1 - theta_v) * (lambda0_final / 2.0) * exp(-lambda0_final * ab),
    1e-300
  )))

  log_lik + log_prior_beta + log_prior_s
}
