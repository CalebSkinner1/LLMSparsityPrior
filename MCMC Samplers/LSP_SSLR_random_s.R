# LSP Spike and Slab Lasso with random sparsity  - MAP Estimation
# adapted from SSLASSO on Cran
# https://cran.r-project.org/web/packages/SSLASSO/index.html
# with data accessed from this github https://github.com/cran/SSLASSO/blob/master/R/SSLASSO.R

# wrapper function -------------------------------------------------------
lsp_random_ssl_map <- function(
  X,
  y,
  weights,
  a_s = 1,
  b_s = NA,
  E_space = c(1, 2, 3),
  c_space = c(0, 1),
  lambda1 = NULL,
  lambda0s = NULL,
  variance = "unknown",
  max_iter = 500,
  eps = 0.001
) {
  n <- nrow(X)
  p <- ncol(X)

  if (is.null(weights)) {
    weights <- rep(1, p)
    E_space <- 0
    c_space <- 0
  }

  if (is.na(b_s)) {
    b_s <- p
  }

  # Standardise X (centre only)
  std_X <- standardize(X)
  XX <- std_X$XX
  mean_y <- mean(y)
  yy <- y - mean_y

  if (is.null(lambda1)) {
    lambda1 <- 1
  }
  if (is.null(lambda0s)) {
    lambda0s <- seq(lambda1, n, length.out = 100)
  }

  # Sigma initialisation
  df <- 3
  sigquant <- 0.9
  sigest <- sd(yy)
  qchi <- qchisq(1 - sigquant, df)
  ncp <- sigest^2 * qchi / df
  min_sigma2 <- sigest^2 / n
  sigma_init <- if (variance == "unknown") sqrt(df * ncp / (df + 2)) else 1.0

  initialbeta <- rep(0, p)
  count_max <- 10L

  # ── Pre-build all v_vec combinations once (memoize eta powers) ──────────
  # For each unique eta, weights^eta and mean(weights^eta) are computed once.
  eta_memo <- list() # keyed by as.character(eta_val)
  all_combos <- vector("list", length(c_space) * length(E_space))
  idx <- 1L
  for (c_val in c_space) {
    for (eta_val in E_space) {
      key <- as.character(eta_val)
      if (is.null(eta_memo[[key]])) {
        w_eta <- weights^eta_val
        eta_memo[[key]] <- list(w_eta = w_eta, mean_w_eta = mean(w_eta))
      }
      memo <- eta_memo[[key]]
      v_vec <- if (c_val == 0) {
        rep(1.0, p)
      } else {
        c_val * (memo$w_eta / memo$mean_w_eta) + (1 - c_val)
      }
      all_combos[[idx]] <- list(c_val = c_val, eta_val = eta_val, v_vec = v_vec)
      idx <- idx + 1L
    }
  }

  map_results <- find_MAP_hyperparams(
    X = XX,
    y = yy,
    weights = weights,
    a_s = a_s,
    b_s = b_s,
    all_combos = all_combos, # pre-built combos passed directly
    lambda1 = lambda1,
    lambda0s = lambda0s,
    initialbeta = initialbeta,
    variance = variance,
    sigma = sigma_init,
    min_sigma2 = min_sigma2,
    eps = eps,
    max_iter = max_iter,
    count_max = count_max
  )

  # Post-processing: recover intercept on original scale
  for (l in seq_along(map_results)) {
    map_results[[l]]$intercept <- mean_y -
      as.numeric(crossprod(std_X$c, map_results[[l]]$beta))
    names(map_results[[l]]$beta) <- colnames(X)
  }

  map_results
}


# ── Helper functions ---------------------------------------------------------

standardize <- function(X) {
  c_vec <- colMeans(X)
  XX <- sweep(X, 2, c_vec, "-")
  list(XX = XX, c = c_vec)
}

expectation_approx <- function(beta_vec, a_s, b_s) {
  p <- length(beta_vec)
  max(0, min(1, (sum(beta_vec != 0) + a_s) / (a_s + b_s + p)))
}

# construct_v_vec is retained for any external callers, but inside the package
# the memoized path in lsp_random_ssl_map() is used instead.
construct_v_vec <- function(weights, c_val, eta) {
  if (c_val == 0) {
    return(rep(1.0, length(weights)))
  }
  w_eta <- weights^eta
  c_val * (w_eta / mean(w_eta)) + (1 - c_val)
}

# Vectorised threshold — unchanged (already fast)
threshold_func_vec <- function(theta_vec, sigma2, lambda1, lambda0, xnorm_vec) {
  p <- length(theta_vec)
  if (lambda1 == lambda0) {
    return(rep(sigma2 * lambda1, p))
  }

  ratio0 <- ((1 - theta_vec) / theta_vec) * (lambda0 / lambda1)
  pstar0 <- 1 / (1 + ratio0)
  lstar0 <- pstar0 * lambda1 + (1 - pstar0) * lambda0
  g0 <- (lstar0 - lambda1)^2 + 2 * xnorm_vec / sigma2 * log(pstar0)

  ifelse(
    g0 > 0,
    sqrt(2 * xnorm_vec * sigma2 * log(1 / pstar0)) + sigma2 * lambda1,
    sigma2 * lstar0
  )
}

update_sigma2 <- function(r) sum(r^2) / length(r)

# Scalar SSL thresholding — unchanged (called per-coordinate, must stay scalar)
SSL_thresholding <- function(
  z,
  beta_val,
  lambda0,
  lambda1,
  theta_j,
  norm,
  delta,
  sigma2
) {
  if (abs(z) <= delta) {
    return(0)
  }
  if (lambda1 == lambda0) {
    lam <- lambda1
  } else {
    ratio <- ((1 - theta_j) / theta_j) *
      (lambda0 / lambda1) *
      exp(-abs(beta_val) * (lambda0 - lambda1))
    pstar_val <- 1 / (1 + ratio)
    lam <- pstar_val * lambda1 + (1 - pstar_val) * lambda0
  }
  temp <- abs(z) - sigma2 * lam
  if (temp > 0) sign(z) * temp / norm else 0
}

# Log-posterior — optimized: fused vectorised prior, fewer allocations
compute_log_posterior <- function(
  y,
  X,
  beta_final,
  sigma_final,
  v_vec,
  lambda1,
  lambda0_final,
  a_s,
  b_s
) {
  n <- length(y)
  p <- length(beta_final)
  sig2 <- sigma_final^2

  s_val <- max(0, min(1, (sum(beta_final != 0) + a_s) / (a_s + b_s + p)))
  # Clamp theta in one pass
  theta_v <- pmin(pmax(s_val * v_vec, 1e-4), 1 - 1e-4)

  # Log-likelihood (residuals computed inline — avoids a named intermediate)
  rss <- sum((y - as.vector(X %*% beta_final))^2)
  log_lik <- -(n / 2) * log(sig2) - rss / (2 * sig2)

  # Log-prior for beta: fused into a single vectorised expression
  ab <- abs(beta_final)
  # pmax on the mixture density guards against log(0); 1e-300 floor is retained
  log_prior_beta <- sum(log(pmax(
    theta_v *
      (lambda1 / 2) *
      exp(-lambda1 * ab) +
      (1 - theta_v) * (lambda0_final / 2) * exp(-lambda0_final * ab),
    1e-300
  )))

  s_safe <- max(min(s_val, 1 - 1e-10), 1e-10)
  log_prior_s <- (a_s - 1) * log(s_safe) + (b_s - 1) * log(1 - s_safe)

  log_lik + log_prior_beta + log_prior_s
}


# ── Core coordinate descent --------------------------------------------------

lsp_ssl_random_descent <- function(
  X,
  y,
  initialbeta,
  variance = "unknown",
  lambda1,
  lambda0s,
  v_vec,
  a_s,
  b_s,
  sigma,
  min_sigma2,
  eps,
  max_iter,
  count_max
) {
  n <- nrow(X)
  p <- ncol(X)
  L <- length(lambda0s)

  xnorm <- colSums(X^2)
  b_mat <- matrix(0, nrow = p, ncol = L)
  sigmas <- rep(NA_real_, L)
  loss <- rep(NA_real_, L)
  iter_vec <- rep(0L, L)
  s_path <- rep(NA_real_, L)

  a_vec <- initialbeta
  newa <- initialbeta
  a_old <- numeric(p)

  e1 <- as.integer(a_vec != 0)
  e2 <- as.integer(a_vec != 0)

  r <- y - as.vector(X %*% a_vec)
  z <- as.vector(crossprod(X, r))

  # Choose computation strategy once
  thres <- min(n, max_iter)
  large_p <- p >= thres # TRUE  → maintain r; FALSE → use XTX cache
  if (!large_p) {
    XTY <- as.vector(crossprod(X, y))
    XTX <- crossprod(X)
  } else {
    XTY <- XTX <- NULL
  }

  sigma2 <- sigma^2
  sigma2_init <- sigma^2
  estimate_sigma <- FALSE
  theta_vec <- pmax(pmin(0.5 * v_vec, 1 - 1e-10), 1e-10)
  delta <- numeric(p)

  for (l in seq_len(L)) {
    lambda0 <- lambda0s[l]

    # Sigma warm-start
    if (l > 1L && variance == "unknown") {
      if (iter_vec[l - 1L] < 100L) {
        estimate_sigma <- TRUE
        sigma2_new <- update_sigma2(r)
        sigma2 <- if (sigma2_new < min_sigma2) {
          estimate_sigma <- FALSE
          sigma2_init
        } else {
          sigma2_new
        }
      } else {
        estimate_sigma <- FALSE
        if (iter_vec[l - 1L] == max_iter) sigma2 <- sigma2_init
      }
    }

    # Theta / delta warm-start
    if (l > 1L) {
      s_val <- expectation_approx(b_mat[, l - 1L], a_s, b_s)
      theta_vec <- pmax(pmin(s_val * v_vec, 1 - 1e-10), 1e-10)
    }

    delta <- threshold_func_vec(theta_vec, sigma2, lambda1, lambda0, xnorm)
    e2 <- pmax(e2, as.integer(abs(z) > delta))

    # active_counter tracks visits to *active* variables only, matching the
    # original intent of "fire every count_max coordinate updates that matter"
    active_counter <- 0L

    # ── Outer while ─────────────────────────────────────────────────────────
    while (iter_vec[l] < max_iter) {
      # ── Active-set pass ──────────────────────────────────────────────────
      while (iter_vec[l] < max_iter) {
        iter_vec[l] <- iter_vec[l] + 1L
        a_old[] <- a_vec

        for (j in seq_len(p)) {
          if (!e1[j]) {
            next
          } # skip inactive variables entirely

          active_counter <- active_counter + 1L

          z[j] <- if (large_p) {
            crossprod(X[, j], r)[[1L]] + xnorm[j] * a_vec[j]
          } else {
            XTY[j] - crossprod(XTX[, j], newa)[[1L]] + xnorm[j] * a_vec[j]
          }

          b_val <- SSL_thresholding(
            z[j],
            a_vec[j],
            lambda0,
            lambda1,
            theta_vec[j],
            xnorm[j],
            delta[j],
            sigma2
          )
          b_mat[j, l] <- b_val

          if (large_p) {
            shift <- b_val - a_vec[j]
            if (shift != 0) r <- r - shift * X[, j]
          } else {
            newa[j] <- b_val
          }

          if (b_val == 0) {
            e1[j] <- 0L
          }

          # Periodic refresh — every count_max *active* updates
          if (active_counter == count_max) {
            if (variance == "unknown" && estimate_sigma) {
              sigma2_new <- update_sigma2(r)
              if (sigma2_new >= min_sigma2) {
                sigma2 <- sigma2_new
              } else {
                sigma2 <- sigma2_init
              }
            }
            s_val <- expectation_approx(b_mat[, l], a_s, b_s)
            theta_vec <- pmax(pmin(s_val * v_vec, 1 - 1e-10), 1e-10)
            delta <- threshold_func_vec(
              theta_vec,
              sigma2,
              lambda1,
              lambda0,
              xnorm
            )
            active_counter <- 0L
          }
        }

        # Sync a_vec
        if (large_p) {
          a_vec <- b_mat[, l]
        } else {
          a_vec <- newa
          newa <- a_vec
        }

        # Convergence check — vectorised denominator with single ifelse
        check_idx <- which(e1 == 1L | a_old != 0)
        converged_active <- TRUE
        if (length(check_idx) > 0L) {
          av <- a_vec[check_idx]
          aov <- a_old[check_idx]
          den <- ifelse(
            av != 0,
            abs(av),
            ifelse(aov != 0, abs(aov), .Machine$double.eps)
          )
          if (any(abs(av - aov) / den > eps)) converged_active <- FALSE
        }
        if (converged_active) break
      }

      # ── Strong-set violation scan ────────────────────────────────────────
      violations_strong <- 0L
      strong_cands <- which(e1 == 0L & e2 == 1L)

      for (j in strong_cands) {
        # Always use residual path; XTX path only where truly beneficial
        z[j] <- if (large_p) {
          crossprod(X[, j], r)[[1L]] + xnorm[j] * a_vec[j]
        } else {
          XTY[j] - crossprod(XTX[, j], a_vec)[[1L]] + xnorm[j] * a_vec[j]
        }

        b_val <- SSL_thresholding(
          z[j],
          a_vec[j],
          lambda0,
          lambda1,
          theta_vec[j],
          xnorm[j],
          delta[j],
          sigma2
        )
        b_mat[j, l] <- b_val

        if (b_val != 0) {
          e1[j] <- 1L
          e2[j] <- 1L
          if (large_p) {
            r <- r - b_val * X[, j]
          }
          a_vec[j] <- b_val
          violations_strong <- violations_strong + 1L

          active_counter <- active_counter + 1L
          if (active_counter == count_max) {
            if (variance == "unknown" && estimate_sigma) {
              sigma2_new <- update_sigma2(r)
              if (sigma2_new >= min_sigma2) {
                sigma2 <- sigma2_new
              } else {
                sigma2 <- sigma2_init
              }
            }
            s_val <- expectation_approx(b_mat[, l], a_s, b_s)
            theta_vec <- pmax(pmin(s_val * v_vec, 1 - 1e-10), 1e-10)
            delta <- threshold_func_vec(
              theta_vec,
              sigma2,
              lambda1,
              lambda0,
              xnorm
            )
            active_counter <- 0L
          }
        }
      }
      if (violations_strong > 0L) {
        next
      }

      # ── Rest violation scan ──────────────────────────────────────────────
      violations_rest <- 0L
      rest_cands <- which(e2 == 0L)

      for (j in rest_cands) {
        z[j] <- if (large_p) {
          crossprod(X[, j], r)[[1L]] + xnorm[j] * a_vec[j]
        } else {
          XTY[j] - crossprod(XTX[, j], a_vec)[[1L]] + xnorm[j] * a_vec[j]
        }

        b_val <- SSL_thresholding(
          z[j],
          a_vec[j],
          lambda0,
          lambda1,
          theta_vec[j],
          xnorm[j],
          delta[j],
          sigma2
        )
        b_mat[j, l] <- b_val

        if (b_val != 0) {
          e1[j] <- 1L
          e2[j] <- 1L
          if (large_p) {
            r <- r - b_val * X[, j]
          }
          a_vec[j] <- b_val
          violations_rest <- violations_rest + 1L

          active_counter <- active_counter + 1L
          if (active_counter == count_max) {
            if (variance == "unknown" && estimate_sigma) {
              sigma2_new <- update_sigma2(r)
              if (sigma2_new >= min_sigma2) {
                sigma2 <- sigma2_new
              } else {
                sigma2 <- sigma2_init
              }
            }
            s_val <- expectation_approx(b_mat[, l], a_s, b_s)
            theta_vec <- pmax(pmin(s_val * v_vec, 1 - 1e-10), 1e-10)
            delta <- threshold_func_vec(
              theta_vec,
              sigma2,
              lambda1,
              lambda0,
              xnorm
            )
            active_counter <- 0L
          }
        }
      }
      if (violations_rest > 0L) {
        next
      }

      # ── Finalise ─────────────────────────────────────────────────────────
      if (!large_p) {
        r <- y - as.vector(X %*% a_vec)
      }
      b_mat[, l] <- a_vec
      loss[l] <- sum(r^2)
      sigmas[l] <- sqrt(sigma2)
      s_path[l] <- expectation_approx(b_mat[, l], a_s, b_s)
      break
    }

    # Fallback if max_iter reached without break
    if (is.na(loss[l])) {
      if (!large_p) {
        r <- y - as.vector(X %*% a_vec)
      }
      b_mat[, l] <- a_vec
      loss[l] <- sum(r^2)
      sigmas[l] <- sqrt(sigma2)
      s_path[l] <- expectation_approx(b_mat[, l], a_s, b_s)
    }
  }

  list(
    beta = b_mat,
    loss = loss,
    iter = iter_vec,
    sigmas = sigmas,
    s_path = s_path
  )
}


# ── Grid-search wrapper ------------------------------------------------------

find_MAP_hyperparams <- function(
  X,
  y,
  weights,
  a_s,
  b_s,
  all_combos, # pre-built list of list(c_val, eta_val, v_vec)
  lambda1,
  lambda0s,
  initialbeta,
  variance,
  sigma,
  min_sigma2,
  eps,
  max_iter,
  count_max
) {
  L <- length(lambda0s)
  n_c <- length(all_combos)

  # Run all (c, eta) paths and store results
  all_paths <- vector("list", n_c)
  for (ci in seq_len(n_c)) {
    combo <- all_combos[[ci]]
    res <- lsp_ssl_random_descent(
      X = X,
      y = y,
      v_vec = combo$v_vec,
      a_s = a_s,
      b_s = b_s,
      variance = variance,
      lambda1 = lambda1,
      lambda0s = lambda0s,
      initialbeta = initialbeta,
      sigma = sigma,
      min_sigma2 = min_sigma2,
      eps = eps,
      max_iter = max_iter,
      count_max = count_max
    )
    all_paths[[ci]] <- list(
      c_val = combo$c_val,
      eta_val = combo$eta_val,
      v_vec = combo$v_vec,
      res = res
    )
  }

  # ── Cross-sectional MAP selection ─────────────────────────────────────────
  # Pre-extract beta/sigma matrices across all paths to avoid repeated indexing
  # inside the double loop.  For each path, cache the L-length sigma vector and
  # the p×L beta matrix so the inner loop only does scoring arithmetic.

  best_path_results <- vector("list", L)
  p <- nrow(all_paths[[1]]$res$beta)

  for (l in seq_len(L)) {
    best_score <- -Inf
    best_params_l <- list(
      lambda0 = lambda0s[l],
      c = NA,
      eta = NA,
      s = NA,
      beta = NULL,
      sigma = NA,
      score = NA
    )

    for (path in all_paths) {
      sigma_l <- path$res$sigmas[l]
      if (is.na(sigma_l) || sigma_l <= 0) {
        next
      }

      beta_l <- path$res$beta[, l]
      if (anyNA(beta_l)) {
        next
      }

      score <- compute_log_posterior(
        y = y,
        X = X,
        beta_final = beta_l,
        sigma_final = sigma_l,
        v_vec = path$v_vec,
        lambda1 = lambda1,
        lambda0_final = lambda0s[l],
        a_s = a_s,
        b_s = b_s
      )

      if (!is.finite(score)) {
        next
      }

      if (score > best_score) {
        best_score <- score
        best_params_l$c <- path$c_val
        best_params_l$eta <- path$eta_val
        best_params_l$s <- path$res$s_path[l]
        best_params_l$beta <- beta_l
        best_params_l$sigma <- sigma_l
        best_params_l$score <- score
      }
    }

    if (is.infinite(best_score)) {
      warning(sprintf("No valid score for lambda0[%d] = %.4f", l, lambda0s[l]))
    }

    best_path_results[[l]] <- best_params_l
  }

  best_path_results
}
