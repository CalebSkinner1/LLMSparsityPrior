# AKI Datasets Server
# In this file, we apply each method on five data sets.

message("loading functions...")
source("aki_analysis_support.R")

# Load Data Sets ---------------------------------------------------------

aki_data_0 <- read_csv(
  "t60_reg_data.csv",
  show_col_types = FALSE
) |>
  select(-pat_id)

# 80 year old
aki_data1 <- aki_data_0 |>
  filter(age > 80) |>
  topK_features(threshold = 0.90)

# female smokers
aki_data2 <- aki_data_0 |>
  filter(current_smoker == 1, gender == 2) |>
  topK_features(threshold = 0.90)

# black men
aki_data3 <- aki_data_0 |>
  filter(race_black == 1, gender == 1) |>
  topK_features(threshold = 0.90)

# liver disease
aki_data4 <- aki_data_0 |>
  filter(liver_dis == 1) |>
  topK_features(threshold = 0.90)

# immunocompromised
aki_data5 <- aki_data_0 |>
  filter(imm_supp == 1) |>
  topK_features(threshold = 0.90)

# Load Weights -----------------------------------------------------------

aki_weights_0 <- read_csv(
  "aki_weights_original_1.csv",
  show_col_types = FALSE
)

aki_weights_continuous <- read_csv(
  "aki_weights_continuous_1.csv",
  show_col_types = FALSE
)

subset_weights <- function(aki_weights_0, aki_data_set) {
  aki_weights_0 |>
    select(value, importance) |>
    filter(value %in% colnames(aki_data_set)) |>
    # ensure in the same order as aki_data columns
    mutate(value = factor(value, levels = colnames(aki_data_set))) |>
    arrange(value)
}

aki_weights1 <- subset_weights(aki_weights_0, aki_data1)
aki_weights2 <- subset_weights(aki_weights_0, aki_data2)
aki_weights3 <- subset_weights(aki_weights_0, aki_data3)
aki_weights4 <- subset_weights(aki_weights_0, aki_data4)
aki_weights5 <- subset_weights(aki_weights_0, aki_data5)

data_weights_list <- list(
  "dataset1" = list("data" = aki_data1, "weights" = aki_weights1),
  "dataset2" = list("data" = aki_data2, "weights" = aki_weights2),
  "dataset3" = list("data" = aki_data3, "weights" = aki_weights3),
  "dataset4" = list("data" = aki_data4, "weights" = aki_weights4),
  "dataset5" = list("data" = aki_data5, "weights" = aki_weights5)
)

aki_cont_weights1 <- subset_weights(aki_weights_continuous, aki_data1)
aki_cont_weights2 <- subset_weights(aki_weights_continuous, aki_data2)
aki_cont_weights3 <- subset_weights(aki_weights_continuous, aki_data3)
aki_cont_weights4 <- subset_weights(aki_weights_continuous, aki_data4)
aki_cont_weights5 <- subset_weights(aki_weights_continuous, aki_data5)

data_weights_continuous_list <- list(
  "dataset1" = list("data" = aki_data1, "weights" = aki_cont_weights1),
  "dataset2" = list("data" = aki_data2, "weights" = aki_cont_weights2),
  "dataset3" = list("data" = aki_data3, "weights" = aki_cont_weights3),
  "dataset4" = list("data" = aki_data4, "weights" = aki_cont_weights4),
  "dataset5" = list("data" = aki_data5, "weights" = aki_cont_weights5)
)

# analysis settings ----------------------------------------------------
tau_range <- c(1, 2, 5, 10)
random_s <- TRUE
fixed_s <- TRUE
sparsity <- 0.01
eta_range <- NULL
iter <- 60000
burn_in <- 10000
outcome <- "creatinine_ratio"
folds <- 5
repetitions <- 10
total_cores <- parallel::detectCores(logical = FALSE)
cores <- min(total_cores, 50)

plan(multicore, workers = cores)
options(future.globals.maxSize = 2000 * 1024^2)

for (j in seq_along(data_weights_list)) {
  message(paste("Begin dataset", j, "..."))

  # Build cross-validation partitions once per dataset
  set.seed(123)
  partitions <- train_test_split(
    data_weights_list[[j]]$data,
    data_weights_list[[j]]$weights,
    n_folds = folds,
    repetitions = repetitions,
    outcome_var = outcome
  )

  set.seed(123)
  cont_partitions <- train_test_split(
    data_weights_continuous_list[[j]]$data,
    data_weights_continuous_list[[j]]$weights,
    n_folds = folds,
    repetitions = repetitions,
    outcome_var = outcome
  )

  # ── Tau-independent models: run ONCE per dataset ──────────────────────
  # message("  running tau-independent models...")
  non_ss_results <- future_map(
    seq_along(partitions),
    ~ {
      train_and_evaluate_non_ss(
        partition = partitions[[.x]],
        cont_partition = cont_partitions[[.x]],
        seed = .x,
        set_eta_range = eta_range,
        random_s = random_s,
        fixed_s = fixed_s
      )
    },
    .options = furrr_options(seed = TRUE)
  ) |>
    transpose() |>
    map(~ .x |> bind_rows())

  # # Tau-dependent models: loop over tau values
  # for (tau in tau_range) {
  #   message(paste("  Begin tau =", tau, "..."))

  #   ss_results <- future_map(
  #     seq_along(partitions),
  #     ~ {
  #       train_and_evaluate_spike_and_slab(
  #         partition = partitions[[.x]],
  #         cont_partition = cont_partitions[[.x]],
  #         seed = .x,
  #         sparsity_type = sparsity_type,
  #         set_tau = tau,
  #         set_eta_range = eta_range,
  #         set_sparsity = sparsity,
  #         set_burn_in = burn_in,
  #         set_iter = iter
  #       )
  #     },
  #     .options = furrr_options(seed = TRUE)
  #   ) |>
  #     transpose() |>
  #     map(~ .x |> bind_rows())

  #   # Write combined results for this (dataset, tau) pair
  #   bind_cols(
  #     non_ss_results$mse,
  #     ss_results$mse |> select(-y)
  #   ) |>
  #     mutate(tau = tau) |>
  #     write_csv(paste0("dataset", j, "tau", tau, "results.csv"))

  #   message(paste("  Completed tau =", tau))
  # }

  non_ss_results$mse |>
    write_csv(paste0("dataset", j, "results.csv"))

  message(paste("Completed dataset", j))
}
