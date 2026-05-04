# AKI Data Application - Five Subsets
#
# Applies all LSP and baseline models to five clinical subgroups derived from
# the AKI dataset. Results are written to one CSV per dataset.

message("loading functions...")
source("AKI Data Application/analysis/aki_analysis_support.R")

# ------------------------------------------------------------------------------
# Load and Filter Data
#
# All subgroups are derived from the same base dataset; near-constant columns
# are removed within each subgroup via topK_features (threshold = 0.90).
# ------------------------------------------------------------------------------

aki_data_0 <- read_csv("t60_reg_data.csv", show_col_types = FALSE) |>
  select(-pat_id)

# Subgroup 1: patients older than 80
aki_data1 <- aki_data_0 |>
  filter(age > 80) |>
  topK_features(threshold = 0.90)
# Subgroup 2: female current smokers
aki_data2 <- aki_data_0 |>
  filter(current_smoker == 1, gender == 2) |>
  topK_features(threshold = 0.90)
# Subgroup 3: Black male patients
aki_data3 <- aki_data_0 |>
  filter(race_black == 1, gender == 1) |>
  topK_features(threshold = 0.90)
# Subgroup 4: patients with liver disease
aki_data4 <- aki_data_0 |>
  filter(liver_dis == 1) |>
  topK_features(threshold = 0.90)
# Subgroup 5: immunocompromised patients
aki_data5 <- aki_data_0 |>
  filter(imm_supp == 1) |>
  topK_features(threshold = 0.90)

# ------------------------------------------------------------------------------
# Load and Align Weights
#
# Two weight sets are used:
#   aki_weights_0             — discretized LLM importance weights (for LSP-SS/SSL)
#   aki_weights_probabilities — probability importance weights (used as direct
#                            prior inclusion probabilities)
#
# subset_weights aligns the weight data frame to the columns present in each
# subgroup dataset after topK_features filtering.
# ------------------------------------------------------------------------------

aki_weights_0 <- read_csv(
  "AKI Data Application/weights/aki_weights_original_1.csv",
  show_col_types = FALSE
)
aki_weights_probabilities <- read_csv(
  "AKI Data Application/weights/aki_weights_probability_1.csv",
  show_col_types = FALSE
)

# Retain and order weights to match the columns of a given subgroup dataset
subset_weights <- function(aki_weights_0, aki_data_set) {
  aki_weights_0 |>
    select(value, importance) |>
    filter(value %in% colnames(aki_data_set)) |>
    mutate(value = factor(value, levels = colnames(aki_data_set))) |>
    arrange(value)
}

aki_weights1 <- subset_weights(aki_weights_0, aki_data1)
aki_weights2 <- subset_weights(aki_weights_0, aki_data2)
aki_weights3 <- subset_weights(aki_weights_0, aki_data3)
aki_weights4 <- subset_weights(aki_weights_0, aki_data4)
aki_weights5 <- subset_weights(aki_weights_0, aki_data5)

aki_prob_weights1 <- subset_weights(aki_weights_probabilities, aki_data1)
aki_prob_weights2 <- subset_weights(aki_weights_probabilities, aki_data2)
aki_prob_weights3 <- subset_weights(aki_weights_probabilities, aki_data3)
aki_prob_weights4 <- subset_weights(aki_weights_probabilities, aki_data4)
aki_prob_weights5 <- subset_weights(aki_weights_probabilities, aki_data5)

# Named lists pairing each subgroup with its respective weight sets
data_weights_list <- list(
  dataset1 = list(data = aki_data1, weights = aki_weights1),
  dataset2 = list(data = aki_data2, weights = aki_weights2),
  dataset3 = list(data = aki_data3, weights = aki_weights3),
  dataset4 = list(data = aki_data4, weights = aki_weights4),
  dataset5 = list(data = aki_data5, weights = aki_weights5)
)

data_weights_probability_list <- list(
  dataset1 = list(data = aki_data1, weights = aki_prob_weights1),
  dataset2 = list(data = aki_data2, weights = aki_prob_weights2),
  dataset3 = list(data = aki_data3, weights = aki_prob_weights3),
  dataset4 = list(data = aki_data4, weights = aki_prob_weights4),
  dataset5 = list(data = aki_data5, weights = aki_prob_weights5)
)

# ------------------------------------------------------------------------------
# Analysis Settings
# ------------------------------------------------------------------------------

outcome <- "creatinine_ratio"
folds <- 5
repetitions <- 10


tau_range <- 1
sparsity <- 0.01
eta_range <- NULL # NULL triggers default prior
iter <- 60000
burn_in <- 10000
random_s <- TRUE
fixed_s <- TRUE

# ------------------------------------------------------------------------------
# Parallel Backend
# ------------------------------------------------------------------------------

total_cores <- parallel::detectCores(logical = FALSE)
cores <- min(total_cores, 50)
plan(multicore, workers = cores)
options(future.globals.maxSize = 2000 * 1024^2)

# ------------------------------------------------------------------------------
# Main Analysis Loop
# ------------------------------------------------------------------------------

for (j in seq_along(data_weights_list)) {
  message(paste("Begin dataset", j, "..."))

  # Build CV partitions once per dataset
  set.seed(123)
  partitions <- train_test_split(
    data_weights_list[[j]]$data,
    data_weights_list[[j]]$weights,
    n_folds = folds,
    repetitions = repetitions,
    outcome_var = outcome
  )

  set.seed(123)
  prob_partitions <- train_test_split(
    data_weights_probability_list[[j]]$data,
    data_weights_probability_list[[j]]$weights,
    n_folds = folds,
    repetitions = repetitions,
    outcome_var = outcome
  )

  non_ss_results <- future_map(
    seq_along(partitions),
    ~ {
      train_and_evaluate_non_ss(
        partition = partitions[[.x]],
        prob_partition = prob_partitions[[.x]],
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

  ss_results <- future_map(
    seq_along(partitions),
    ~ {
      train_and_evaluate_spike_and_slab(
        partition = partitions[[.x]],
        prob_partition = prob_partitions[[.x]],
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

  bind_cols(
    non_ss_results$mse,
    ss_results$mse |> select(-y)
  ) |>
    mutate(tau = tau) |>
    write_csv(paste0("dataset", j, "results.csv"))

  message("Completed dataset ", j)
}
