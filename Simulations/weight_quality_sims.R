# R script for running simulations in LLM Sparsity Prior for Robust Feature Selection
# source("Simulations/weight_quality_support.R")
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
eta_range <- seq(by = 1, from = 1, to = 10)
confidence_range <- c(0, 1)
iter <- 30000
burn_in <- 5000
random_s <- TRUE
fixed_s <- FALSE

# select grid of weight agreements to evaluate
phi_range <- c(0.5, 0.6, 0.7, 0.75, seq(0.8, 1.0, by = 0.01))

total_cores <- parallel::detectCores(logical = FALSE)
cores <- min(total_cores, 50)

plan(multisession, workers = cores)
options(future.globals.maxSize = 2000 * 1024^2)

# run simulations ----------------------------------------------------------
for (n in n_range) {
  message(paste0("Generating data and baseline models for n = ", n, "..."))

  cached_baselines <- future_map(
    1:n_replications,
    function(seed_idx) {
      baseline_data_sim_function(seed = seed_idx, n = n)
    },
    .options = furrr_options(seed = TRUE)
  )

  for (phi in phi_range) {
    # select weight quality phi
    message(paste0("   evaluating weights for phi = ", phi))
    weights <- generate_weights(phi, true_gamma, categories = 5)
    l1_agreement <- l1_weight_agreement(true_gamma, weights)
    l2_agreement <- l2_weight_agreement(true_gamma, weights)
    pairwise_agreement <- pairwise_weight_agreement(true_gamma, weights)
    roc_agreement <- ROC_weight_agreement(true_gamma, weights)

    file_name <- paste0("weights", phi, "_n", n, ".csv")

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
