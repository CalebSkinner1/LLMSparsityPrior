# LLM Sparsity Prior (this time place a prior on the sparsity)
# Spike-and-Slab Lasso for Regression with a unique constant prior inclusion probability for each covariate

# compute model prior conditional on s
compute_log_prior_gamma <- function(gamma, s, u) {
  theta <- pmax(1e-10, pmin(s * u, 1 - 1e-10)) # compute theta vector

  # compute log probability of gamma conditional on s
  log_prob <- sum(gamma * log(theta) + (1 - gamma) * log(1 - theta))

  log_prob
}

# run spike and slab lasso gibbs sampler for regression.
lsp_random_ssl_gibbs_sampler <- function(
  X,
  y,
  weights = NULL,
  c = NA,
  eta = 0,
  a_sigma,
  b_sigma,
  lambda_0 = 50,
  lambda_1 = 1,
  a_s = 1,
  b_s = NA,
  s_proposal_sigma = 1,
  iter = 10000,
  burn_in = 5000,
  thin = 1,
  init_weights = TRUE,
  return_samples = TRUE
) {
  if (is.null(weights)) {
    c <- 0 # if no weights, then clearly there is no confidence in them
    eta <- 0
  }

  p <- ncol(X)
  n <- nrow(X)

  # number of values of eta and c in grid
  K <- length(c) * length(eta)

  # all c, eta values
  cross_eta_c <- expand.grid(eta = eta, c = c)

  # create space for u_mat (u*sparsity = theta)
  u_mat <- matrix(0, nrow = nrow(cross_eta_c), ncol = p)

  # eta and c are fixed
  if (length(c) == 1 & length(eta) == 1) {
    if(eta == 0 | c == 0){ 
      init_weights <- FALSE
      u_mat[1,] <- rep(1, p)
    } else {
      # create vector for prior model probability
      u_mat[1,] <- c * (weights^eta) / mean(weights^eta) + (1 - c)
    }
  } else { # eta and c are random
    for(k in 1:nrow(cross_eta_c)){
      c_k <- cross_eta_c$c[k]
      eta_k <- cross_eta_c$eta[k]

      u_mat[k, ] <- c_k * (weights^eta_k) / mean((weights^eta_k)) + (1 - c_k)
    }
  }

  # per recommendation of rockova-george
  if (is.na(b_s)) {
    b_s <- p
  }

  # number of models left after thinning/burn_in
  n_keep <- ceiling((iter - burn_in) / thin)

  if(return_samples){
    # create space for gamma, beta, invsigma^2, etc. (acc_store removed for gamma)
    gam_store <- matrix(0, nrow = n_keep, ncol = p)
    beta_store <- matrix(0, nrow = n_keep, ncol = p + 1)
    invsigma_2_store <- rep(0, n_keep)
    eta_store <- rep(0, n_keep)
    c_store <- rep(0, n_keep)
    s_store <- rep(0, n_keep)
    acc_s_store <- rep(0, n_keep)
  } else {
    # reduction for memory: only store means
    gam_mean <- rep(0, p)
    beta_mean <- rep(0, p + 1)
    invsigma_2_mean <- 0
    eta_mean <- 0
    c_mean <- 0
    s_mean <- 0
    acc_s_mean <- 0
  }

  # generate initial values of gamma, s, sigma^2 (in a smart way)
  gam_current <- rep(0, p)
  s_current <- a_s / (a_s + b_s)

  if (init_weights) {
    # initialize by taking largest p*prior_mean_s weights from LLM
    gam_current[order(weights, decreasing = TRUE)[
      1:max(2, ceiling(s_current * p))
    ]] <- 1
  } else {
    # initialize by taking largest p*prior_mean_s correlations with y
    gam_current[order(abs(cor(X, as.vector(y))), decreasing = TRUE)[
      1:max(2, ceiling(s_current * p))
    ]] <- 1
  }

  # Initialize beta using ridge
  fit <- glmnet::glmnet(X, y, alpha = 0, lambda = 1)
  beta_current <- as.vector(coef(fit))
  rm(fit)

  invsigma_2_current <- as.numeric(1 / var(as.vector(y)))
  v_current <- rep(1, p) 
  X_aug <- cbind(1, X)
  
  # find current index of eta and c in discrete uniform grid
  if(K > 1){
    c_eta_idx_current <- sample(K, 1)
  } else {
    c_eta_idx_current <- 1
  }
  
  # begin iterations
  for (i in 1:iter) {
    
    sigma_current <- as.numeric(sqrt(1 / invsigma_2_current))
    theta_current <- pmax(1e-10, pmin(s_current * u_mat[c_eta_idx_current, ], 1 - 1e-10))

    # 1. Update Gamma 
    beta_no_int <- beta_current[-1]
    
    log_p1 <- log(theta_current) + log(lambda_1) - (lambda_1 * abs(beta_no_int) / sigma_current)
    log_p0 <- log(1 - theta_current) + log(lambda_0) - (lambda_0 * abs(beta_no_int) / sigma_current)
    
    prob_gam <- 1 / (1 + exp(log_p0 - log_p1))
    gam_current <- rbinom(p, 1, prob_gam)

    # 2. Update Latent Scales (v) using Inverse Gaussian
    lambda_current <- ifelse(gam_current == 1, lambda_1, lambda_0)
    mu_inv_v <- sqrt( (lambda_current^2 * sigma_current^2) / (beta_no_int^2 + 1e-10) )
    
    inv_v_current <- statmod::rinvgauss(p, mean = mu_inv_v, shape = lambda_current^2)
    v_current <- 1 / inv_v_current

    # 3. Update Beta
    D_inv <- diag(c(0, inv_v_current)) 
    
    Q <- crossprod(X_aug) + D_inv
    chol_Q <- chol(Q)
    invQ <- chol2inv(chol_Q)
    
    mu_beta <- as.vector(invQ %*% crossprod(X_aug, y)) 
    beta_current <- MASS::mvrnorm(1, mu = mu_beta, Sigma = invQ / as.numeric(invsigma_2_current))

    # 4. Update invsigma_2
    shape_sigma <- (n + p) / 2 + a_sigma
    rate_sigma <- b_sigma + 
                  0.5 * sum((as.vector(y) - as.vector(X_aug %*% beta_current))^2) + 
                  0.5 * sum((beta_current[-1]^2) * inv_v_current)
    
    invsigma_2_current <- rgamma(1, shape = shape_sigma, rate = rate_sigma)

    # 5. draw new s via random walk on logit scale
    logit_s <- log(s_current / (1 - s_current))
    logit_s_new <- rnorm(1, mean = logit_s, sd = s_proposal_sigma)
    s_new <- 1 / (1 + exp(-logit_s_new))

    if (max(s_new * u_mat[c_eta_idx_current,]) >= 1) {
      accept_s <- FALSE
    } else {
      # compute posterior ratio for s
      log_prior_s_new <- dbeta(s_new, a_s, b_s, log = TRUE)
      log_prior_s_old <- dbeta(s_current, a_s, b_s, log = TRUE)

      log_lik_s_new <- compute_log_prior_gamma(gam_current, s_new, u_mat[c_eta_idx_current,])
      log_lik_s_old <- compute_log_prior_gamma(gam_current, s_current, u_mat[c_eta_idx_current,])

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

    # 6. draw new c and eta based on gamma
    # unnormalized log probability for all K states of eta and c
    W <- numeric(K)
    for(k in 1:K){
      theta_k <- pmin(pmax(s_current * u_mat[k,], 1e-4), 1 - 1e-4)
      W[k] <- sum(gam_current * log(theta_k) + (1 - gam_current) * log(1 - theta_k))
    }
    # normalize probabilities
    # log-sum-exp trick to prevent NaN underflow
    pi_eta_c <- exp(W - max(W)) / sum(exp(W - max(W)))
    c_eta_idx_current <- sample(K, 1, prob = pi_eta_c)

    # store parameters
    if (i > burn_in && (i - burn_in) %% thin == 0) {
      c_val <- cross_eta_c$c[c_eta_idx_current]
      if(return_samples){
        store_i <- (i - burn_in) / thin # index

        gam_store[store_i, ] <- gam_current
        beta_store[store_i, ] <- beta_current
        invsigma_2_store[store_i] <- invsigma_2_current
        c_store[store_i] <- c_val
        eta_store[store_i] <- if (c_val == 0) 0 else cross_eta_c$eta[c_eta_idx_current]
        s_store[store_i] <- s_current
        acc_s_store[store_i] <- accept_s
      } else {
        gam_mean <- gam_mean + gam_current / n_keep
        beta_mean <- beta_mean + beta_current / n_keep
        invsigma_2_mean <- invsigma_2_mean + invsigma_2_current / n_keep
        c_mean <- c_mean + c_val / n_keep
        eta_mean <- if(c_val == 0) eta_mean else eta_mean + cross_eta_c$eta[c_eta_idx_current] / n_keep
        s_mean <- s_mean + s_current / n_keep
        acc_s_mean <- acc_s_mean + accept_s / n_keep
      }
    }
  }

  if(return_samples){
    list(
      "beta" = beta_store,
      "gamma" = gam_store,
      "invsigma_2" = invsigma_2_store,
      "eta" = eta_store,
      "c" = c_store,
      "s" = s_store,
      "acc_s" = acc_s_store
    )
  } else {
    list(
      "beta" = beta_mean,
      "gamma" = gam_mean,
      "invsigma_2" = invsigma_2_mean,
      "eta" = eta_mean,
      "c" = c_mean,
      "s" = s_mean,
      "acc_s" = acc_s_mean
    )
  }
}