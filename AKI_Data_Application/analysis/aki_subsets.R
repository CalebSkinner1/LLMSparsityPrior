# AKI Data Application - Five Subsets
#
# Applies all LSP and baseline models to five clinical subgroups derived from
# the AKI dataset. Results are written to one CSV per dataset.

message("loading functions...")
source("AKI_Data_Application/analysis/aki_analysis_support.R")

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
# Two weight sets are used for each subgroup:
#   - discretized LLM importance weights (for LSP-SS/SSL)
#   - probability importance weights (used as direct prior inclusion
#     probabilities)
#
# subset_weights aligns the weight data frame to the columns present in each
# subgroup dataset after topK_features filtering.
# ------------------------------------------------------------------------------

aki_weights_1 <- read_csv(
  "AKI_Data_Application/weights/aki_weights_original_1.csv",
  show_col_types = FALSE
)
aki_weights_probabilities_1 <- read_csv(
  "AKI_Data_Application/weights/aki_weights_probabilities_1.csv",
  show_col_types = FALSE
)

aki_weights_2 <- read_csv(
  "AKI_Data_Application/weights/aki_weights_original_female_smoker.csv",
  show_col_types = FALSE
)
aki_weights_probabilities_2 <- read_csv(
  "AKI_Data_Application/weights/aki_weights_probabilities_female_smoker.csv",
  show_col_types = FALSE
)

aki_weights_3 <- read_csv(
  "AKI_Data_Application/weights/aki_weights_original_black_men.csv",
  show_col_types = FALSE
)
aki_weights_probabilities_3 <- read_csv(
  "AKI_Data_Application/weights/aki_weights_probabilities_black_men.csv",
  show_col_types = FALSE
)

aki_weights_4 <- read_csv(
  "AKI_Data_Application/weights/aki_weights_original_liver_disease.csv",
  show_col_types = FALSE
)
aki_weights_probabilities_4 <- read_csv(
  "AKI_Data_Application/weights/aki_weights_probabilities_liver_disease.csv",
  show_col_types = FALSE
)

aki_weights_5 <- read_csv(
  "AKI_Data_Application/weights/aki_weights_original_immunocompromised.csv",
  show_col_types = FALSE
)
aki_weights_probabilities_5 <- read_csv(
  "AKI_Data_Application/weights/aki_weights_probabilities_immunocompromised.csv",
  show_col_types = FALSE
)

aki_weights1 <- subset_weights(aki_weights_1, aki_data1)
aki_weights2 <- subset_weights(aki_weights_2, aki_data2)
aki_weights3 <- subset_weights(aki_weights_3, aki_data3)
aki_weights4 <- subset_weights(aki_weights_4, aki_data4)
aki_weights5 <- subset_weights(aki_weights_5, aki_data5)

aki_prob_weights1 <- subset_weights(aki_weights_probabilities_1, aki_data1)
aki_prob_weights2 <- subset_weights(aki_weights_probabilities_2, aki_data2)
aki_prob_weights3 <- subset_weights(aki_weights_probabilities_3, aki_data3)
aki_prob_weights4 <- subset_weights(aki_weights_probabilities_4, aki_data4)
aki_prob_weights5 <- subset_weights(aki_weights_probabilities_5, aki_data5)

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

tau <- 1
sparsity <- 0.01
eta_range <- NULL # NULL triggers default prior
iter <- 60000
burn_in <- 10000
thin <- 5
random_s <- TRUE
fixed_s <- FALSE

# ------------------------------------------------------------------------------
# Parallel Backend
# ------------------------------------------------------------------------------

total_cores <- parallel::detectCores(logical = FALSE)
cores <- min(total_cores, 12)
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

  # Same partitions with the naive LLM inclusion probabilities swapped in
  prob_partitions <- map(partitions, function(partition) {
    partition$weights <- data_weights_probability_list[[j]]$weights
    partition
  })

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
        random_s = random_s,
        fixed_s = fixed_s,
        set_tau = tau,
        set_eta_range = eta_range,
        set_sparsity = sparsity,
        set_burn_in = burn_in,
        set_iter = iter,
        set_thin = thin,
        return_coverage = TRUE
      )
    },
    .options = furrr_options(seed = TRUE)
  ) |>
    transpose()

  ss_mse <- ss_results$mse |> bind_rows()
  ss_coverage <- ss_results$coverage |>
    bind_rows(.id = "partition") |>
    mutate(
      partition = as.integer(partition),
      rep = (partition - 1L) %/% folds + 1L,
      fold = (partition - 1L) %% folds + 1L
    )

  # per-replication coverage
  ss_coverage |>
    group_by(method, rep) |>
    summarize(
      n_intervals = n(),
      coverage = mean(coverage),
      mean_width = mean(width),
      median_width = median(width),
      .groups = "drop"
    ) |>
    mutate(dataset = j, tau = tau, nominal = 0.95, .before = 1) |>
    write_csv(paste0("dataset", j, "coverage_by_rep.csv"))

  bind_cols(
    non_ss_results$mse,
    ss_mse |> select(-y)
  ) |>
    mutate(tau = tau) |>
    write_csv(paste0("dataset", j, "results.csv"))

  message("Completed dataset ", j)
}
