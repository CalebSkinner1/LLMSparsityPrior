# R script for using server to run simulations

source("server_reg_support.R")
# source("Simulations/Regression/server - use for multiple cores/server_reg_support.R")

# simulation settings -----------------------------------------------------

p <- 1000
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
tau <- 2
sparsity <- 0.01
random_s <- FALSE
fixed_s <- TRUE
# select grid of weight agreements to evaluate
phi_range <- c(0.5, 0.6, 0.7, 0.75, seq(0.8, 1.0, by = 0.01))

total_cores <- parallel::detectCores(logical = FALSE)
cores <- min(total_cores, 50)
weight_confidence <- 0.5

plan(multisession, workers = cores)
options(future.globals.maxSize = 2000 * 1024^2)

# run simulations ----------------------------------------------------------
for (phi in phi_range) {
  message(str_c("starting simulations for phi = ", phi))
  weights <- generate_weights(phi, true_gamma)
  l1_agreement <- l1_weight_agreement(true_gamma, weights)
  l2_agreement <- l2_weight_agreement(true_gamma, weights)
  pairwise_agreement <- pairwise_weight_agreement(true_gamma, weights)

  for (n in c(250, 500)) {
    message(str_c("running n = ", n, "..."))
    file_name <- paste0("weights", phi, "_n", n, ".csv")

    sim_results <- future_map(
      1:n_replications,
      function(seed_idx) {
        sim_function(seed = seed_idx, n = n, weights = weights)
      },
      .options = furrr_options(seed = TRUE)
    ) |>
      map_dfr(~ .x |> bind_rows(.id = "method"), .id = "sim_id") |>
      mutate(
        l1_agreement = l1_agreement,
        l2_agreement = l2_agreement,
        pairwise_agreement = pairwise_agreement,
        n = n
      )

    # write result
    write_csv(sim_results, file_name)
  }
}
