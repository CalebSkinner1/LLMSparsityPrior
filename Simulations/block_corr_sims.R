# Simulation Driver — Weight Quality Study (Block Correlation)
#
# Evaluates LSP model performance across a grid of sample sizes (n_range) and
# synthetic weight quality levels (phi_range) under a block-diagonal covariance
# structure: covariates are partitioned into equal-sized blocks with
# exchangeable correlation Xcorr within each block and zero correlation
# between blocks. Signal positions are randomly permuted in each replicate,
# with the weight vector permuted identically, so results do not depend on
# how the signals line up with the blocks. For each n, baseline models (no LLM weights)
# are fit once and cached; LSP and LLM-Lasso models are then fit for
# each phi level reusing the same datasets. Results are written to one
# CSV per (phi, n) combination.
#
# Depends on: Simulations/weight_quality_support.R

source("Simulations/weight_quality_support.R")

# ------------------------------------------------------------------------------
# Simulation Settings
# ------------------------------------------------------------------------------

# Data Generating Process
p <- 1000
n_range <- c(100, 250)
s <- 20
# Canonical ordering (signals last). Weights are generated against this
# ordering and permuted per replicate inside sim_function.
true_gamma <- c(rep(0, p - s), rep(1, s))
randomize_beta <- TRUE
effect_size <- 1
Xvar <- 1
Xcorr <- 0.7
block_size <- 50
y_sd <- 1
cov_mat <- block_cor_mat(p, block_size = block_size, rho = Xcorr)

# Sampler hyperparameters
a_sigma <- 1
b_sigma <- 1
tau <- 1
sparsity <- 0.01
eta_range <- seq(from = 0, to = 10, by = 1)
iter <- 30000
burn_in <- 5000

# Model variants to run
random_s <- TRUE
fixed_s <- FALSE

# Simulation grid
phi_range <- seq(0.5, 1.0, by = 0.1)
n_replications <- 500

# ------------------------------------------------------------------------------
# Parallel Backend
# ------------------------------------------------------------------------------

total_cores <- parallel::detectCores(logical = FALSE)
cores <- min(total_cores, 25)
plan(multicore, workers = cores)
options(future.globals.maxSize = 2000 * 1024^2)

# ------------------------------------------------------------------------------
# Simulation Loop
# ------------------------------------------------------------------------------

for (n in n_range) {
  # Generate all replicate datasets and fit baseline models once,
  # then reuse across phi levels.
  message(paste0("Generating data and baseline models for n = ", n, "..."))

  cached_baselines <- future_map(
    1:n_replications,
    function(seed_idx) {
      baseline_data_sim_function(
        seed = seed_idx,
        n = n,
        randomize_beta = randomize_beta
      )
    },
    .options = furrr_options(seed = TRUE)
  )

  for (phi in phi_range) {
    # Construct a synthetic weight vector achieving L1 agreement level phi
    message(paste0("   evaluating weights for phi = ", phi))
    weights <- generate_weights(phi, true_gamma, categories = 5)
    l1_agreement <- l1_weight_agreement(true_gamma, weights)
    l2_agreement <- l2_weight_agreement(true_gamma, weights)
    pairwise_agreement <- pairwise_weight_agreement(true_gamma, weights)
    roc_agreement <- ROC_weight_agreement(true_gamma, weights)

    file_name <- paste0("block", Xcorr, "_weights", phi, "_n", n, ".csv")

    # Fit LSP and LLM-Lasso models across all replicates in parallel;
    # bind results into a long tibble indexed by simulation ID and method.
    sim_results <- future_map(
      cached_baselines,
      function(baseline) {
        sim_function(baseline_fits = baseline, weights = weights)
      },
      .options = furrr_options(seed = TRUE)
    ) |>
      map_dfr(~ .x |> bind_rows(.id = "method"), .id = "sim_id") |>
      mutate(
        l1_agreement = l1_agreement,
        l2_agreement = l2_agreement,
        pairwise_agreement = pairwise_agreement,
        roc_agreement = roc_agreement,
        n = n
      )
    # write result
    write_csv(sim_results, file_name)
  }
}
