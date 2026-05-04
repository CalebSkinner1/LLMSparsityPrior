# script that runs low data regimes on aki analysis on server
message("loading functions...")
source("aki_analysis_support.R")

# load data --------------------------------------------------------------

aki_data <- read_csv(
  "t60_reg_data.csv",
  show_col_types = FALSE
) |>
  select(-pat_id) |>
  # subset on age to reduce data size
  filter(age > 80) |> # 80 year old
  # select features with at least a 90% uniqueness rate
  topK_features(threshold = 0.90)

# load weights -----------------------------------------------------------

aki_weights <- read_csv(
  "aki_weights_original_1.csv",
  show_col_types = FALSE
) |>
  select(value, importance) |>
  filter(value %in% colnames(aki_data)) |>
  # ensure in the same order as aki_data columns
  mutate(value = factor(value, levels = colnames(aki_data))) |>
  arrange(value)

aki_weights_continuous <- read_csv(
  "aki_weights_continuous_1.csv",
  show_col_types = FALSE
) |>
  select(value, importance) |>
  filter(value %in% colnames(aki_data)) |>
  # ensure in the same order as aki_data columns
  mutate(value = factor(value, levels = colnames(aki_data))) |>
  arrange(value)

# analysis settings ----------------------------------------------------
tau <- 1
sparsity <- 0.01
sparsity_type <- "random"
eta_range <- NULL
iter <- 60000
burn_in <- 10000
outcome <- "creatinine_ratio"
folds <- 5
repetitions <- 10
n_range <- c(50, 100, 150, 200)
total_cores <- parallel::detectCores(logical = FALSE)
cores <- min(total_cores, 50)

plan(multicore, workers = cores)
options(future.globals.maxSize = 2000 * 1024^2)

message("running low data analysis...")
for (n in n_range) {
  set.seed(123) # ensure partitions are the same each time
  # partition the data into 10 5-fold cross validation for evaluation
  partitions <- train_test_split(
    aki_data,
    aki_weights,
    n_folds = folds,
    repetitions = repetitions,
    outcome_var = outcome,
    n = n
  )

  message("    running baseline models...")
  baseline_results <- future_map(
    seq_along(partitions),
    ~ {
      train_and_evaluate_baselines(
        partition = partitions[[.x]],
        seed = .x,
        sparsity_type = sparsity_type,
        set_tau = tau,
        set_sparsity = sparsity,
        set_burn_in = burn_in,
        set_iter = iter
      )
    },
    .options = furrr_options(seed = TRUE)
  ) |>
    transpose() |>
    map(~ .x |> bind_rows())

  message("    running llm methods...")
  llm_methods_results <- future_map(
    seq_along(partitions),
    ~ {
      train_and_evaluate_random_eta(
        partition = partitions[[.x]],
        seed = .x,
        sparsity_type = sparsity_type,
        set_tau = tau,
        set_eta_range = eta_range,
        set_sparsity = sparsity,
        set_burn_in = burn_in,
        set_iter = iter
      )
    },
    .options = furrr_options(seed = TRUE)
  ) |>
    transpose() |>
    map(~ .x |> bind_rows())

  message("    running continuous weights...")
  set.seed(123)
  cont_partitions <- train_test_split(
    aki_data,
    aki_weights_continuous,
    n_folds = folds,
    repetitions = repetitions,
    outcome_var = outcome,
    n = n
  )

  continuous_results <- future_map(
    seq_along(cont_partitions),
    ~ {
      train_and_evaluate_continuous_weights(
        partition = cont_partitions[[.x]],
        seed = .x,
        set_tau = tau,
        set_burn_in = burn_in,
        set_iter = iter
      )
    },
    .options = furrr_options(seed = TRUE)
  ) |>
    transpose() |>
    map(~ .x |> bind_rows())

  csv_file_name <- paste0("low_data_results_n", n, ".csv")

  # write mse
  bind_cols(
    baseline_results$mse,
    llm_methods_results$mse |> select(-y),
    continuous_results$mse |> select(-y)
  ) |>
    write_csv(csv_file_name)
  message(paste("   Successfully saved:", csv_file_name))
}
