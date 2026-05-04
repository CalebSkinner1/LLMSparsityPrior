# AKI Analysis - Server
# script that runs aki analysis on server
# R script for running AKI Application in LLM Sparsity Prior for Robust Feature Selection
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

read_weights_file <- function(file_path) {
  read_csv(
    paste0(file_path, ".csv"),
    show_col_types = FALSE
  ) |>
    select(value, importance) |>
    filter(value %in% colnames(aki_data)) |>
    # ensure in the same order as aki_data columns
    mutate(value = factor(value, levels = colnames(aki_data))) |>
    arrange(value)
}

weight_names <- c(
  str_c("aki_weights_original_", seq(1, 5, by = 1)),
  str_c("aki_weights_task_", seq(1, 5, by = 1)),
  str_c("aki_weights_reasoning_", seq(1, 5, by = 1)),
  str_c("aki_weights_collinearity_", seq(1, 5, by = 1)),
  str_c("aki_weights_rubric10_", seq(1, 5, by = 1))
)

aki_weights_list <- map(weight_names, read_weights_file)
names(aki_weights_list) <- weight_names

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
total_cores <- parallel::detectCores(logical = FALSE)
cores <- min(total_cores, 50)

plan(multicore, workers = cores)
options(future.globals.maxSize = 2000 * 1024^2)

message("running weight sensitivity...")
for (weight_name in names(aki_weights_list)) {
  set.seed(123) # ensure partitions are the same each time
  # partition the data into 10 5-fold cross validation for evaluation
  partitions <- train_test_split(
    aki_data,
    aki_weights_list[[weight_name]],
    n_folds = folds,
    repetitions = repetitions,
    outcome_var = outcome
  )

  weight_sensitivity_results <- future_map(
    seq_along(partitions),
    ~ {
      train_and_evaluate_random_eta(
        partition = partitions[[.x]],
        seed = .x,
        sparsity_type = sparsity_type,
        set_tau = tau,
        set_sparsity = sparsity,
        set_eta_range = eta_range,
        set_confidence_range = confidence_range,
        set_burn_in = burn_in,
        set_iter = iter
      )
    },
    .options = furrr_options(seed = TRUE)
  ) |>
    transpose() |>
    map(~ .x |> bind_rows())

  csv_file_name <- paste0("weight_sensitivity_results_", weight_name, ".csv")

  # write mse
  write_csv(weight_sensitivity_results$mse, csv_file_name)
  message(paste("   Successfully saved:", csv_file_name))
}
