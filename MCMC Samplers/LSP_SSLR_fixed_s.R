# LSP Spike and Slab Lasso with Fixed sparsity  - MAP Estimation
# adapted from SSLASSO on Cran
# https://cran.r-project.org/web/packages/SSLASSO/index.html
# with data accessed from this github https://github.com/cran/SSLASSO/blob/master/R/SSLASSO.R

# wrapper function -------------------------------------------------------
lsp_fixed_ssl_map <- function(
  X,
  y,
  weights = NULL,
  s_fixed = 0.05,
  E_space = NULL,
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
  } else if (is.null(E_space)) {
    eta_max <- 0
    step_size <- 1
    # Calculate initial bound to ensure it starts < 1
    theta_bound <- s_fixed * max(weights)^eta_max / mean(weights^eta_max)

    while (eta_max <= 20) {
      eta_max <- eta_max + step_size # step forward
      theta_bound <- s_fixed * max(weights)^eta_max / mean(weights^eta_max)

      # if threshold is crossed, backtrack with smaller steps
      if (theta_bound >= 1) {
        # Step back to the last safe value
        eta_max <- eta_max - step_size

        # Decrease the step size for finer searching
        if (step_size == 1) {
          step_size <- 0.1
        } else if (step_size == 0.1) {
          step_size <- 0.01
        } else {
          break
        }
      }
    }
    # generate eta space
    E_space <- seq(0, eta_max, length.out = 11)
    rm(eta_max)
  }

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

  # ── Pre-build all (c, eta) v_vecs once; memoize eta powers ──────────────
  eta_memo <- list()
  all_combos <- vector("list", length(E_space))
  idx <- 1L
  for (eta_val in E_space) {
    key <- as.character(eta_val)
    if (is.null(eta_memo[[key]])) {
      w_eta <- weights^eta_val
      eta_memo[[key]] <- list(w_eta = w_eta, mean_w_eta = mean(w_eta))
    }
    memo <- eta_memo[[key]]
    v_vec <- if (eta_val == 0) {
      rep(1.0, p)
    } else {
      (memo$w_eta / memo$mean_w_eta)
    }
    all_combos[[idx]] <- list(eta_val = eta_val, v_vec = v_vec)
    idx <- idx + 1L
  }

  map_results <- find_MAP_fixed_hyperparams(
    X = XX,
    y = yy,
    weights = weights,
    s_fixed = s_fixed,
    all_combos = all_combos,
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

construct_v_vec <- function(weights, eta) {
  if (eta == 0) {
    return(rep(1.0, length(weights)))
  }
  w_eta <- weights^eta
  (w_eta / mean(w_eta))
}

threshold_func_vec <- function(theta_vec, sigma2, lambda1, lambda0, xnorm_vec) {
  if (lambda1 == lambda0) {
    return(rep(sigma2 * lambda1, length(theta_vec)))
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

update_sigma2 <- function(r) sum(r^2) / (length(r) + 2)

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

# Fixed-s log-posterior: s_fixed is a constant so the s prior is omitted.
# term1 + term2 fused into one vectorised expression to cut allocations.
compute_log_posterior_fixed <- function(
  y,
  X,
  beta_final,
  sigma_final,
  v_vec,
  lambda1,
  lambda0_final,
  s_fixed
) {
  n <- length(y)
  sig2 <- sigma_final^2

  # theta_vec is constant for the entire run
  if (length(s_fixed) == p) {
    theta_v <- pmax(pmin(s_fixed, 1 - 1e-10), 1e-10)
  } else {
    theta_v <- pmax(pmin(s_fixed * v_vec, 1 - 1e-10), 1e-10)
  }

  rss <- sum((y - as.vector(X %*% beta_final))^2)
  log_lik <- -(n / 2) * log(sig2) - rss / (2 * sig2)

  ab <- abs(beta_final)
  log_prior_beta <- sum(log(pmax(
    theta_v *
      (lambda1 / 2) *
      exp(-lambda1 * ab) +
      (1 - theta_v) * (lambda0_final / 2) * exp(-lambda0_final * ab),
    1e-300
  )))

  log_lik + log_prior_beta
}


# ── Core coordinate descent --------------------------------------------------

lsp_ssl_fixed_descent <- function(
  X,
  y,
  initialbeta,
  variance = "unknown",
  lambda1,
  lambda0s,
  v_vec,
  s_fixed,
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

  a_vec <- initialbeta
  newa <- initialbeta
  a_old <- numeric(p)

  e1 <- as.integer(a_vec != 0)
  e2 <- as.integer(a_vec != 0)

  r <- y - as.vector(X %*% a_vec)
  z <- as.vector(crossprod(X, r))

  # Computation strategy decided once
  thres <- min(n, max_iter)
  large_p <- p >= thres
  if (!large_p) {
    XTY <- as.vector(crossprod(X, y))
    XTX <- crossprod(X)
  } else {
    XTY <- XTX <- NULL
  }

  sigma2 <- sigma^2
  sigma2_init <- sigma^2
  estimate_sigma <- FALSE

  # theta_vec is constant for the entire run — computed once here
  if (length(s_fixed) == p) {
    theta_vec <- pmax(pmin(s_fixed, 1 - 1e-10), 1e-10)
  } else {
    theta_vec <- pmax(pmin(s_fixed * v_vec, 1 - 1e-10), 1e-10)
  }

  delta <- numeric(p)

  for (l in seq_len(L)) {
    lambda0 <- lambda0s[l]

    # Sigma warm-start
    if (l > 1L && variance == "unknown") {
      if (iter_vec[l - 1L] < 100L) {
        estimate_sigma <- TRUE
        sigma2_new <- update_sigma2(r)
        if (sigma2_new < min_sigma2) {
          sigma2 <- sigma2_init
          estimate_sigma <- FALSE
        } else {
          sigma2 <- sigma2_new
        }
      } else {
        estimate_sigma <- FALSE
        if (iter_vec[l - 1L] == max_iter) sigma2 <- sigma2_init
      }
    }

    # delta must update on every new lambda0 (lambda0 always changes)
    delta <- threshold_func_vec(theta_vec, sigma2, lambda1, lambda0, xnorm)
    e2 <- pmax(e2, as.integer(abs(z) > delta))

    # active_counter: only increments on active-variable visits
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
          }

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

          # Periodic refresh every count_max *active* updates.
          # theta_vec is constant, so delta is only recomputed when sigma2
          # actually changes — saves a full length-p vectorised call otherwise.
          if (active_counter == count_max) {
            if (variance == "unknown" && estimate_sigma) {
              sigma2_new <- update_sigma2(r)
              if (sigma2_new >= min_sigma2) {
                sigma2 <- sigma2_new
                delta <- threshold_func_vec(
                  theta_vec,
                  sigma2,
                  lambda1,
                  lambda0,
                  xnorm
                )
              } else {
                sigma2 <- sigma2_init
              }
            }
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

        # Convergence check — single-pass vectorised denominator
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
                delta <- threshold_func_vec(
                  theta_vec,
                  sigma2,
                  lambda1,
                  lambda0,
                  xnorm
                )
              } else {
                sigma2 <- sigma2_init
              }
            }
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
                delta <- threshold_func_vec(
                  theta_vec,
                  sigma2,
                  lambda1,
                  lambda0,
                  xnorm
                )
              } else {
                sigma2 <- sigma2_init
              }
            }
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
      break
    }

    # Fallback
    if (is.na(loss[l])) {
      if (!large_p) {
        r <- y - as.vector(X %*% a_vec)
      }
      b_mat[, l] <- a_vec
      loss[l] <- sum(r^2)
      sigmas[l] <- sqrt(sigma2)
    }
  }

  list(beta = b_mat, loss = loss, iter = iter_vec, sigmas = sigmas)
}


# ── Grid-search wrapper ------------------------------------------------------

find_MAP_fixed_hyperparams <- function(
  X,
  y,
  weights,
  s_fixed,
  all_combos,
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

  # Run all (c, eta) paths
  all_paths <- vector("list", n_c)
  for (ci in seq_len(n_c)) {
    combo <- all_combos[[ci]]
    res <- lsp_ssl_fixed_descent(
      X = X,
      y = y,
      v_vec = combo$v_vec,
      s_fixed = s_fixed,
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
      eta_val = combo$eta_val,
      v_vec = combo$v_vec,
      res = res
    )
  }

  # Cross-sectional MAP selection
  best_path_results <- vector("list", L)

  for (l in seq_len(L)) {
    best_score <- -Inf
    best_params_l <- list(
      lambda0 = lambda0s[l],
      eta = NA,
      s = s_fixed,
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

      score <- compute_log_posterior_fixed(
        y = y,
        X = X,
        beta_final = beta_l,
        sigma_final = sigma_l,
        v_vec = path$v_vec,
        lambda1 = lambda1,
        lambda0_final = lambda0s[l],
        s_fixed = s_fixed
      )

      # zero-inflate eta == 0, to give baseline model 50% prior probability
      if (path$eta_val != 0) {
        score <- score - log(length(all_paths) - 1)
      }

      if (!is.finite(score)) {
        next
      }

      if (score > best_score) {
        best_score <- score
        best_params_l$eta <- path$eta_val
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
