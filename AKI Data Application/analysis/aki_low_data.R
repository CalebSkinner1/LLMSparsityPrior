# AKI Data Application — Low-Data Regime Analysis
#
# Fits baseline and LSP models on the over-80 subgroup across a range of
# training set sizes (n_range), using repeated stratified cross-validation.
# LSP and naive weight variants are evaluated in parallel.
# Results are written to one CSV per training size n.

message("loading functions...")
source("AKI Data Application/analysis/aki_analysis_support.R")

# ------------------------------------------------------------------------------
# Load and Filter Data
#
# Analysis is restricted to patients older than 80. Near-constant columns are
# removed via topK_features (threshold = 0.90).
# ------------------------------------------------------------------------------

aki_data <- read_csv("t60_reg_data.csv", show_col_types = FALSE) |>
  select(-pat_id) |>
  filter(age > 80) |>
  topK_features(threshold = 0.90)

# ------------------------------------------------------------------------------
# Load and Align Weights
#
# Weights are filtered to the columns retained after topK_features and ordered
# to match aki_data's column layout.
# ------------------------------------------------------------------------------

aki_weights <- read_csv(
  "AKI Data Application/weights/aki_weights_original_1.csv",
  show_col_types = FALSE
) |>
  select(value, importance) |>
  filter(value %in% colnames(aki_data)) |>
  mutate(value = factor(value, levels = colnames(aki_data))) |>
  arrange(value)

aki_weights_probabilities <- read_csv(
  "AKI Data Application/weights/aki_weights_probabilities_1.csv",
  show_col_types = FALSE
) |>
  select(value, importance) |>
  filter(value %in% colnames(aki_data)) |>
  mutate(value = factor(value, levels = colnames(aki_data))) |>
  arrange(value)


# ------------------------------------------------------------------------------
# Analysis Settings
# ------------------------------------------------------------------------------

outcome <- "creatinine_ratio"
folds <- 5
repetitions <- 10
n_range <- c(100, 150, 200) # training set sizes to evaluate

tau <- 1
sparsity <- 0.01
eta_range <- NULL # NULL triggers automatic grid search in samplers
iter <- 60000
burn_in <- 10000
fixed_s <- FALSE
random_s <- TRUE

# ------------------------------------------------------------------------------
# Parallel Backend
# ------------------------------------------------------------------------------
total_cores <- parallel::detectCores(logical = FALSE)
cores <- min(total_cores, 50)
plan(multicore, workers = cores)
options(future.globals.maxSize = 2000 * 1024^2)

# ------------------------------------------------------------------------------
# Low-Data Regime Loop
# ------------------------------------------------------------------------------

message("running low data analysis...")
for (n in n_range) {
  message("  n = ", n)
  # Partitions for LSP and naive weights use the same seed so
  # that train/test splits are aligned when results are combined
  set.seed(123)
  partitions <- train_test_split(
    aki_data,
    aki_weights,
    n_folds = folds,
    repetitions = repetitions,
    outcome_var = outcome,
    n = n
  )

  set.seed(123)
  prob_partitions <- train_test_split(
    aki_data,
    aki_weights_probabilities,
    n_folds = folds,
    repetitions = repetitions,
    outcome_var = outcome,
    n = n
  )

  message("    running baseline models...")
  baseline_results <- future_map(
    seq_along(partitions),
    ~ train_and_evaluate_baselines(
      partition = partitions[[.x]],
      seed = .x,
      fixed_s = fixed_s,
      random_s = random_s,
      set_tau = tau,
      set_sparsity = sparsity,
      set_burn_in = burn_in,
      set_iter = iter
    ),
    .options = furrr_options(seed = TRUE)
  ) |>
    transpose() |>
    map(~ bind_rows(.x))

  message("    running llm methods...")
  llm_methods_results <- future_map(
    seq_along(partitions),
    ~ train_and_evaluate_random_eta(
      partition = partitions[[.x]],
      seed = .x,
      fixed_s = fixed_s,
      random_s = random_s,
      set_tau = tau,
      set_eta_range = eta_range,
      set_sparsity = sparsity,
      set_burn_in = burn_in,
      set_iter = iter
    ),
    .options = furrr_options(seed = TRUE)
  ) |>
    transpose() |>
    map(~ bind_rows(.x))

  message("    running naive approach (probabilities)...")
  naive_results <- future_map(
    seq_along(prob_partitions),
    ~ train_and_evaluate_probability_weights(
      partition = prob_partitions[[.x]],
      seed = .x,
      set_tau = tau,
      set_burn_in = burn_in,
      set_iter = iter
    ),
    .options = furrr_options(seed = TRUE)
  ) |>
    transpose() |>
    map(~ bind_rows(.x))

  csv_file_name <- paste0("low_data_results_n", n, ".csv")

  # write mse
  bind_cols(
    baseline_results$mse,
    llm_methods_results$mse |> select(-y),
    naive_results$mse |> select(-y)
  ) |>
    write_csv(csv_file_name)

  message("  Completed n = ", n)
}
