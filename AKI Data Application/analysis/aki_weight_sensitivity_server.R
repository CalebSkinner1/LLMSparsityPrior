# AKI Data Application — Weight Sensitivity Analysis
#
# Evaluates LSP model performance across 25 LLM weight sets spanning five
# prompt variants (original, task, reasoning, collinearity, rubric10), each
# with five independent LLM runs. The over-80 subgroup is used throughout.
# Results are written to one CSV per weight set.
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
# Load Weight Sets
#
# Each weight file is filtered to the columns retained after topK_features and
# ordered to match aki_data's column layout. All 25 files are loaded into a
# named list for iteration.
# ------------------------------------------------------------------------------

read_weights_file <- function(file_path) {
  read_csv(
    paste0(file_path, ".csv"),
    show_col_types = FALSE
  ) |>
    select(value, importance) |>
    filter(value %in% colnames(aki_data)) |>
    mutate(value = factor(value, levels = colnames(aki_data))) |>
    arrange(value)
}

weight_names <- c(
  str_c("aki_weights_original_", 1:5),
  str_c("aki_weights_task_", 1:5),
  str_c("aki_weights_reasoning_", 1:5),
  str_c("aki_weights_collinearity_", 1:5),
  str_c("aki_weights_rubric10_", 1:5)
)

aki_weights_list <- map(
  str_c("AKI Data Application/analysis/weights/", weight_names),
  read_weights_file
)
names(aki_weights_list) <- weight_names

# ------------------------------------------------------------------------------
# Analysis Settings
# ------------------------------------------------------------------------------

outcome <- "creatinine_ratio"
folds <- 5
repetitions <- 10

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
# Weight Sensitivity Loop
# ------------------------------------------------------------------------------

message("Running weight sensitivity analysis...")
for (weight_name in names(aki_weights_list)) {
  message("  Weight set: ", weight_name)

  set.seed(123)
  partitions <- train_test_split(
    aki_data,
    aki_weights_list[[weight_name]],
    n_folds = folds,
    repetitions = repetitions,
    outcome_var = outcome
  )

  future_map(
    seq_along(partitions),
    ~ train_and_evaluate_random_eta(
      partition = partitions[[.x]],
      seed = .x,
      fixed_s = fixed_s,
      random_s = random_s,
      set_tau = tau,
      set_sparsity = sparsity,
      set_eta_range = eta_range,
      set_burn_in = burn_in,
      set_iter = iter
    ),
    .options = furrr_options(seed = TRUE)
  ) |>
    transpose() |>
    map(~ bind_rows(.x)) |>
    pluck("mse") |>
    write_csv(paste0("weight_sensitivity_results_", weight_name, ".csv"))

  message("  Completed: ", weight_name)
}
