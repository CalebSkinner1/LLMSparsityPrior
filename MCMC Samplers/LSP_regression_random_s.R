# LLM Sparsity Prior (this time place a prior on the sparsity)
# Discrete Spike and Slab for Regression with a unique constant prior inclusion probability for each covariate

# compute model prior conditional on s
compute_log_prior_gamma <- function(gamma, s, u) {
  t <- s * u # compute t_j vector

  # compute log probability of gamma conditional on s
  log_prob <- sum(gamma * log(t) + (1 - gamma) * log(1 - t))

  log_prob
}

# compute unnormalized log-posterior density (marginalizing out beta and sigma)
# (Z is (1, X)), does not include some constants that will be cancelled in log_acceptance_rate
# rswr (random sparsity weights regression)
rswr_log_posterior <- function(
  Z,
  Z_gram,
  y,
  tau,
  gamma,
  a_sigma,
  b_sigma,
  s,
  u,
  n
) {
  n_gam <- ncol(Z) # cardinality of selected design matrix

  Q <- Z_gram + diag(n_gam) / tau
  cholQ <- chol(Q)
  log_detQ <- 2 * sum(log(diag(cholQ)))

  model_prior <- compute_log_prior_gamma(gamma, s, u)

  n_gam /
    2 *
    log(tau) -
    .5 * log_detQ -
    (n / 2 + a_sigma) *
      log(
        sum(y^2) - crossprod(y, Z %*% chol2inv(cholQ) %*% t(Z)) %*% y + b_sigma
      ) +
    model_prior
}

# compute log acceptance rate: log(p(gamma_new|data)) - log(p(gamma_old|data))
rswr_log_acceptance_rate <- function(
  Z_old,
  Z_old_gram,
  Z_new,
  Z_new_gram,
  y,
  gamma_new,
  gamma_old,
  tau,
  a_sigma,
  b_sigma,
  s,
  u,
  n
) {
  rswr_log_posterior(
    Z_new,
    Z_new_gram,
    y,
    tau,
    gamma_new,
    a_sigma,
    b_sigma,
    s,
    u,
    n
  ) -
    rswr_log_posterior(
      Z_old,
      Z_old_gram,
      y,
      tau,
      gamma_old,
      a_sigma,
      b_sigma,
      s,
      u,
      n
    )
}

rswr_gibbs_sampler <- function(
  X,
  y,
  weights = NULL,
  c = NA,
  a_sigma,
  b_sigma,
  tau,
  a_s = 1,
  b_s = NA,
  iter = 10000,
  burn_in = 5000,
  thin = 1,
  prob_add = 1 / 3,
  prob_delete = 1 / 3,
  init_weights = TRUE
) {
  if (is.null(weights)) {
    c <- 0 # if no weights, then clearly there is no confidence in them
  }

  p <- ncol(X)
  n <- nrow(X)

  if (c == 0) {
    u <- rep(1, p)
    init_weights <- FALSE
  } else {
    # create vector for prior model probability
    u <- c * weights / mean(weights) + (1 - c)
  }

  # per recommendation of rockova-george
  if (is.na(b_s)) {
    b_s <- p / 10
  }

  # number of models left after thinning/burn_in
  n_keep <- ceiling((iter - burn_in) / thin)

  # create space for gamma, beta, invsigma^2, acc
  gam_list <- matrix(0, nrow = n_keep, ncol = p)
  beta_list <- matrix(0, nrow = n_keep, ncol = p + 1)
  invsigma_2_list <- rep(0, n_keep)
  acc_list <- rep(0, n_keep)

  # generate initial values of gamma, s, sigma^2 (in a smart way)
  gam_current <- rep(0, p)

  s_current <- a_s / (a_s + b_s)

  if (init_weights) {
    # initialize by taking largest p*prior_mean_s weights from LLM
    gam_current[order(weights, decreasing = TRUE)[
      1:max(2, ceiling(s_current * p))
    ]] <- 1
  } else {
    # initialize by taking largest p*prior_mean_s correlations with y (marginal correlation screening)
    gam_current[order(abs(cor(X, y)), decreasing = TRUE)[
      1:max(2, ceiling(s_current * p))
    ]] <- 1
  }

  fit <- glmnet::glmnet(
    X[, which(gam_current == 1)],
    y,
    alpha = 0,
    lambda = 1
  )
  beta_current <- as.vector(coef(fit))

  invsigma_2_current <- 1 /
    (1 /
      n *
      sum((y - cbind(1, X[, which(gam_current == 1)]) %*% beta_current)^2))

  # begin iterations
  for (i in 1:iter) {
    # propose a candidate gamma
    gam_prop <- gam_current
    selected_gam <- which(gam_prop == 1)
    removed_gam <- which(gam_prop == 0)

    Z_old <- cbind(1, X[, selected_gam])
    Z_old_gram <- crossprod(Z_old)

    # employ ADS random search to edit gam_prop
    unif_gam <- runif(1)
    if (length(selected_gam) == 0) {
      unif_gam <- 0.5
    } # ensure that if none are selected, we will add

    log_prop_ratio <- 0
    current_model_size <- sum(gam_prop)

    # with prob_delete, randomly remove one gamma
    if (unif_gam < prob_delete || length(removed_gam) == 0) {
      chosen <- sample(selected_gam)[1]
      gam_prop[chosen] <- 0
      log_prop_ratio <- log(current_model_size) -
        log(p - current_model_size + 1)
    } else if (unif_gam < prob_delete + prob_add) {
      # with prob_add, randomly add one gamma
      chosen <- sample(removed_gam)[1]
      gam_prop[chosen] <- 1
      log_prop_ratio <- log(p - current_model_size) -
        log(current_model_size + 1)
    } else {
      # else swap
      chosen1 <- sample(removed_gam)[1]
      chosen2 <- sample(selected_gam)[1]
      gam_prop[chosen1] <- 1
      gam_prop[chosen2] <- 0
    }
    selected_gam <- which(gam_prop == 1)

    Z_new <- cbind(1, X[, selected_gam])
    Z_new_gram <- crossprod(Z_new)

    # compute log acceptance rate
    logacc <- rswr_log_acceptance_rate(
      Z_old,
      Z_old_gram,
      Z_new,
      Z_new_gram,
      y,
      gam_prop,
      gam_current,
      tau,
      a_sigma,
      b_sigma,
      s_current,
      u,
      n
    ) +
      log_prop_ratio
    if (log(runif(1)) < logacc) {
      # insert gamma draw
      gam_current <- gam_prop

      # draw beta vector
      chol_mat <- chol(
        invsigma_2_current *
          Z_new_gram +
          (invsigma_2_current / tau) * diag(ncol(Z_new_gram))
      )
      invQ <- chol2inv(chol_mat)
      l <- invsigma_2_current * crossprod(Z_new, y)
      beta_gamma <- MASS::mvrnorm(1, invQ %*% l, invQ)
      beta_current <- rep(0, p + 1)
      beta_current[c(1, selected_gam + 1)] <- beta_gamma

      # draw invsigma_2
      invsigma_2_current <- rgamma(
        1,
        shape = n / 2 + a_sigma,
        rate = 1 / 2 * sum((y - Z_new %*% beta_gamma)^2) + b_sigma
      )

      # count acceptances
      acc <- 1
    } else {
      # gam_current is maintained

      # draw beta vector
      chol_mat <- chol(
        invsigma_2_current *
          Z_old_gram +
          (invsigma_2_current / tau) * diag(ncol(Z_old_gram))
      )
      invQ <- chol2inv(chol_mat)
      l <- invsigma_2_current * crossprod(Z_old, y)
      beta_gamma <- MASS::mvrnorm(1, invQ %*% l, invQ)
      beta_current <- rep(0, p + 1)
      beta_current[c(1, which(gam_current == 1) + 1)] <- beta_gamma

      # draw invsigma_2
      invsigma_2_current <- rgamma(
        1,
        shape = n / 2 + a_sigma,
        rate = 1 / 2 * sum((y - Z_old %*% beta_gamma)^2) + b_sigma
      )

      # count acceptances
      acc <- 0
    }

    # draw new s via random walk on logit scale
    logit_s <- log(s_current / (1 - s_current))
    proposal_sigma <- 0.5
    logit_s_new <- rnorm(1, mean = logit_s, sd = proposal_sigma)
    s_new <- 1 / (1 + exp(-logit_s_new))

    if (max(s_new * u) >= 1) {
      accept_s <- FALSE
    } else {
      # compute posterior ratio for s
      log_prior_s_new <- dbeta(s_new, a_s, b_s, log = TRUE)
      log_prior_s_old <- dbeta(s_current, a_s, b_s, log = TRUE)

      log_lik_s_new <- compute_log_prior_gamma(gam_current, s_new, u)
      log_lik_s_old <- compute_log_prior_gamma(gam_current, s_current, u)

      # jacobian adjustment
      log_jacobian_new <- log(s_new) + log(1 - s_new)
      log_jacobian_old <- log(s_current) + log(1 - s_current)

      # total Metrpolis-Hastings ratio for s proposal
      log_ratio_s <- (log_lik_s_new + log_prior_s_new + log_jacobian_new) -
        (log_lik_s_old + log_prior_s_old + log_jacobian_old)

      if (log(runif(1)) < log_ratio_s) {
        s_current <- s_new
        accept_s <- TRUE
      } else {
        accept_s <- FALSE
      }
    }

    # store parameters
    if (i > burn_in && (i - burn_in) %% thin == 0) {
      store_i <- (i - burn_in) / thin # index

      gam_list[store_i, ] <- gam_current
      beta_list[store_i, ] <- beta_current
      invsigma_2_list[store_i] <- invsigma_2_current
      acc_list[store_i] <- acc
    }
  }

  list(
    "beta" = beta_list,
    "gamma" = gam_list,
    "invsigma_2" = invsigma_2_list,
    "accs" = acc_list
  )
}
