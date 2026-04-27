# LSP Spike and Slab Lasso with random sparsity  - MAP Estimation
# adapted from SSLASSO on CRAN
# https://cran.r-project.org/web/packages/SSLASSO/index.html
# with source code accessed from this github https://github.com/cran/SSLASSO/blob/master/R/SSLASSO.R

dll_path <- file.path(getwd(), paste0("LSP_SSL/lsp_ssl", .Platform$dynlib.ext))
if ("lsp_ssl" %in% names(getLoadedDLLs())) {
  dyn.unload(dll_path)
}
dyn.load(dll_path)

stopifnot(is.loaded("SSL_gaussian", PACKAGE = "lsp_ssl"))

# LSP for Spike-and-Slab LASSO
lsp_ssl_map <- function(
  X,
  y,
  weights = NULL,
  E_space = NULL,
  eta_zero_mass = 0.5, # prior mass on eta = 0
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
  if (any(is.na(y)) | any(is.na(X))) {
    stop("Missing data (NA's) detected. Eliminate missing data before calling.")
  }

  XX <- scale(X, center = TRUE, scale = FALSE)
  p <- ncol(XX)
  yy <- y - mean(y)
  n <- length(yy)

  # validate sparsity: scalar or p-dimensional
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

  # validate weights; default to uniform
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

  # build default E_space
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
    # validate given E_space
    E_space <- as.numeric(E_space)
    if (length(E_space) < 1) {
      stop("E_space must contain at least one eta value")
    }
  }

  n_eta <- length(E_space)

  if (missing(lambda0)) {
    lambda0 <- seq(1, n, length = nlambda)
    lambda1 <- lambda0[1]
  } else {
    nlambda <- length(lambda0)
    if (missing(lambda1)) lambda1 <- lambda0[1]
  }

  monotone <- sum((lambda0[-1] - lambda0[-nlambda]) > 0)
  if (monotone != nlambda - 1) {
    stop("lambda0 must be a monotone increasing sequence")
  }
  if (lambda1 > min(lambda0)) {
    stop("lambda1 must be smaller than lambda0")
  }
  if (missing(b_s)) {
    b_s <- p
  }

  df <- 3
  sigquant <- 0.9
  sigest <- sd(yy)
  qchi <- qchisq(1 - sigquant, df)
  ncp <- sigest^2 * qchi / df
  min_sigma2 <- sigest^2 / n
  if (variance == "unknown" && missing(sigma)) {
    sigma <- sqrt(df * ncp / (df + 2))
  }

  # helper function that computes v_vec from weights and a given eta
  make_v_vec <- function(eta) {
    w_eta <- weights^eta
    w_eta / mean(w_eta)
  }

  # run coordinate descent once per eta value
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

  # For each lambda0, select the eta that maximizes the log-posterior
  bb <- matrix(0.0, p, nlambda)
  iter <- integer(nlambda)
  thetas <- numeric(nlambda)
  sigmas <- numeric(nlambda)
  best_eta <- numeric(nlambda) # stores the selected eta value per lambda0

  for (l in seq_len(nlambda)) {
    # Use the baseline sigma as a common reference for all eta comparisons
    sigma_ref <- eta_results[[which(E_space == 0)]]$sigmas[l]
    # Fall back to first run if 0 is not in E_space
    if (length(sigma_ref) == 0 || is.nan(sigma_ref)) {
      sigma_ref <- eta_results[[1]]$sigmas[l]
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

        # add zero-inflated discrete uniform prior on eta
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

    # catch to ensure that one eta is always selected (default is eta == 0 if fails to converge)
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

  if (warn & any(iter == max.iter)) {
    warning(
      "Algorithm did not converge for lambda0: ",
      paste(lambda0[iter == max.iter], collapse = ", ")
    )
  }
  if (iter[nlambda] == max.iter) {
    warning("Algorithm did not converge at the last lambda0 value.")
  }

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
      penalty = penalty,
      lambda1 = lambda1,
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

# compute_log_posterior --------------------------------------------------

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

  ab <- abs(beta_final)
  log_prior_beta <- sum(log(pmax(
    theta_v *
      (lambda1 / 2.0) *
      exp(-lambda1 * ab) +
      (1 - theta_v) * (lambda0_final / 2.0) * exp(-lambda0_final * ab),
    1e-300
  )))

  s_safe <- max(min(s_val, 1.0 - 1e-10), 1e-10)
  log_prior_s <- (a_s - 1.0) * log(s_safe) + (b_s - 1.0) * log(1.0 - s_safe)

  log_lik + log_prior_beta + log_prior_s
}
