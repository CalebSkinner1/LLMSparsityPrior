# AKI Eta Sensitivity Analysis

message("loading functions...")
source("aki_analysis_support.R")

# load data --------------------------------------------------------------

aki_data <- read_csv(
  "t60_reg_data.csv",
  show_col_types = FALSE
) |>
  select(-pat_id) |>
  filter(age > 80) |> # 80 year old
  # select features with at least a 90% uniqueness rate
  topK_features(threshold = 0.90)

aki_weights <- read_csv(
  "aki_weights_original_1.csv",
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
eta_range <- seq(by = 1, from = 1, to = 12)
confidence_range <- 1
iter <- 60000
burn_in <- 10000
outcome <- "creatinine_ratio"
folds <- 5
repetitions <- 10
total_cores <- parallel::detectCores(logical = FALSE)
cores <- min(total_cores, 50)

# partition the data into 10 5-fold cross validation for evaluation
set.seed(123)
partitions <- train_test_split(
  aki_data,
  aki_weights,
  n_folds = folds,
  repetitions = repetitions,
  outcome_var = outcome
)

plan(multicore, workers = cores)
options(future.globals.maxSize = 2000 * 1024^2)

# run hyperparameter sweep for LSP and LLM-Lasso
for (eta in eta_range) {
  # select separation parameter eta
  message(str_c("starting simulations for eta = ", eta))

  # select confidence c
  for (c in confidence_range) {
    message(str_c("running c = ", c, "..."))
    file_name <- paste0("eta", eta, "_c", c, ".csv")

    aki_reg_results <- future_map(
      seq_along(partitions),
      ~ {
        train_and_evaluate_fixed_eta(
          partition = partitions[[.x]],
          seed = .x,
          sparsity_type = sparsity_type,
          set_tau = tau,
          set_sparsity = sparsity,
          set_eta = eta,
          set_confidence = c,
          set_burn_in = burn_in,
          set_iter = iter
        )
      },
      .options = furrr_options(seed = TRUE)
    ) |>
      transpose() |>
      map(~ .x |> bind_rows())

    # write mse
    write_csv(aki_reg_results$mse, file_name)
  }
}
