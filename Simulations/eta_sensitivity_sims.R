# Rscript for running simulations in LLM Sparsity Prior for Robust Feature Selection

source("weight_quality_support.R")

# simulation settings -----------------------------------------------------

p <- 1000
n_range <- c(100, 250)
s <- 20
true_gamma <- c(rep(0, p - s), rep(1, s))
effect_size <- 1
n_replications <- 500
Xvar <- 1
Xcorr <- 0.5
y_sd <- 1
cov_mat <- simstudy::genCorMat(p, cors = rep(Xcorr, choose(p, 2)))
a_sigma <- 1
b_sigma <- 1
tau <- 1
sparsity <- 0.01
eta_range <- seq(by = 1, from = 1, to = 20)
iter <- 30000
burn_in <- 5000
random_s <- TRUE

# select grid of weight agreements to evaluate
phi_range <- c(0.7, 0.8, 0.9)

total_cores <- parallel::detectCores(logical = FALSE)
cores <- min(total_cores, 50)

plan(multisession, workers = cores)
options(future.globals.maxSize = 2000 * 1024^2)

# run simulations ----------------------------------------------------------
for (phi in phi_range) {
  # select weight quality phi
  message(paste0("   evaluating weights for phi = ", phi))
  weights <- generate_weights(phi, true_gamma, categories = 5)
  l1_agreement <- l1_weight_agreement(true_gamma, weights)

  for (eta in eta_range) {
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
    # write result
    write_csv(eta_sensitivity_results, file_name)
  }
}
