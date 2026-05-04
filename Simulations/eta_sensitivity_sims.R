# Simulation Driver — Eta Sensitivity Study
#
# Evaluates LSP model performance across a grid of weight quality levels
# (phi_range) and fixed eta values (eta_range), holding sample size constant.
# For each phi, a synthetic weight vector is constructed; models are then fit
# at each eta value across all replicates in parallel. Results are written to
# one CSV per (phi, eta) combination.
#
# Depends on: Simulations/weight_quality_support.R
source("Simulations/weight_quality_support.R")

# ------------------------------------------------------------------------------
# Simulation Settings
# ------------------------------------------------------------------------------

# Data Generating Process
p <- 1000
n <- 100
s <- 20
true_gamma <- c(rep(0, p - s), rep(1, s))
effect_size <- 1
Xvar <- 1
Xcorr <- 0.5
y_sd <- 1
cov_mat <- simstudy::genCorMat(p, cors = rep(Xcorr, choose(p, 2)))

# Sampler hyperparameters
a_sigma <- 1
b_sigma <- 1
tau <- 1
sparsity <- 0.01
iter <- 30000
burn_in <- 5000

# Model variants to run
random_s <- TRUE
fixed_s <- FALSE

# Simulation grid
phi_range <- c(0.8, 0.9)
eta_range <- seq(1, 20, by = 1)
n_replications <- 500


# ------------------------------------------------------------------------------
# Parallel Backend
# ------------------------------------------------------------------------------

total_cores <- parallel::detectCores(logical = FALSE)
cores <- min(total_cores, 50)
plan(multicore, workers = cores)
options(future.globals.maxSize = 2000 * 1024^2)

# ------------------------------------------------------------------------------
# Simulation Loop
# ------------------------------------------------------------------------------

for (phi in phi_range) {
  message(paste0("Evaluating phi = ", phi))

  weights <- generate_weights(phi, true_gamma, categories = 5)
  l1_agreement <- l1_weight_agreement(true_gamma, weights)

  for (eta in eta_range) {
    message("    Evaluating eta = ", eta)

    file_name <- paste0("weights", phi, "eta_", eta, ".csv")

    eta_sensitivity_results <- future_map(
      1:n_replications,
      function(seed) {
        eta_sensitivity_function(seed = seed, weights = weights, set_eta = eta)
      },
      .options = furrr_options(seed = TRUE)
    ) |>
      map_dfr(~ .x |> bind_rows(.id = "method"), .id = "sim_id") |>
      mutate(
        l1_agreement = l1_agreement,
        n = n
      )
    write_csv(eta_sensitivity_results, file_name)
  }
}
