# AKI Data Application — Low-Data Regime Analysis
#
# Fits baseline and LSP models on each of five clinical subgroups across a
# range of training set sizes (n_range), using repeated stratified cross-validation.
# LSP and naive weight variants are evaluated in parallel.
# Results are written to one CSV per dataset, n combination

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
# Weights are filtered to the columns retained after topK_features and ordered
# to match each aki_data's column layout.
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
n_range <- c(100, 125, 150, 200) # training set sizes to evaluate

tau <- 1
sparsity <- 0.01
eta_range <- NULL # NULL triggers automatic grid search in samplers
iter <- 60000
burn_in <- 10000
thin <- 5
fixed_s <- FALSE
random_s <- TRUE
compute_coverage <- TRUE
coverage_level <- 0.95

# ------------------------------------------------------------------------------
# Parallel Backend
# ------------------------------------------------------------------------------
total_cores <- parallel::detectCores(logical = FALSE)
cores <- min(total_cores, 12)
plan(multicore, workers = cores)
options(future.globals.maxSize = 2000 * 1024^2)

# ------------------------------------------------------------------------------
# Feasible Training Sizes
#
# With v-fold CV the training set holds roughly (folds - 1) / folds of the
# subgroup. train_test_split subsamples when n < nrow(train_data).
# Other combinations are skipped instead.
# ------------------------------------------------------------------------------

max_train_size <- function(data, n_folds) {
  floor(nrow(data) * (n_folds - 1) / n_folds)
}

message("Feasible training sizes by subgroup:")
for (j in seq_along(data_weights_list)) {
  dat <- data_weights_list[[j]]$data
  message(
    "  dataset ",
    j,
    ": ",
    nrow(dat),
    " obs, max training size ",
    max_train_size(dat, folds),
    " -> n in {",
    paste(n_range[n_range <= max_train_size(dat, folds)], collapse = ", "),
    "}"
  )
}

# ------------------------------------------------------------------------------
# Low-Data Regime Loop
# ------------------------------------------------------------------------------

message("running low data analysis...")


for (j in seq_along(data_weights_list)) {
  message("Begin dataset ", j, " ...")

  aki_data <- data_weights_list[[j]]$data
  aki_weights <- data_weights_list[[j]]$weights
  aki_prob_weights <- data_weights_probability_list[[j]]$weights

  max_n <- max_train_size(aki_data, folds)

  for (n in n_range) {
    if (n > max_n) {
      message(
        "  skipping n = ",
        n,
        " (dataset ",
        j,
        " training set holds at most ",
        max_n,
        ")"
      )
      next
    }

    message("  n = ", n)
    # Partitions are built once; the naive weights are swapped in so that
    # train/test splits are aligned when results are combined
    set.seed(123)
    partitions <- train_test_split(
      aki_data,
      aki_weights,
      n_folds = folds,
      repetitions = repetitions,
      outcome_var = outcome,
      n = n
    )

    prob_partitions <- map(partitions, function(partition) {
      partition$weights <- aki_prob_weights
      partition
    })

    # Realized training size, for the record: equals n whenever subsampling
    # actually occurred in every fold.
    realized_n <- min(map_int(partitions, ~ nrow(.x$X_train_scaled)))
    if (realized_n != n) {
      message(
        "    note: realized training size is ",
        realized_n,
        " (requested ",
        n,
        ")"
      )
    }

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
        set_iter = iter,
        set_thin = thin,
        return_coverage = compute_coverage,
        set_level = coverage_level
      ),
      .options = furrr_options(seed = TRUE)
    ) |>
      transpose()

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
        set_iter = iter,
        set_thin = thin,
        return_coverage = compute_coverage,
        set_level = coverage_level
      ),
      .options = furrr_options(seed = TRUE)
    ) |>
      transpose()

    message("    running naive approach (probabilities)...")
    naive_results <- future_map(
      seq_along(prob_partitions),
      ~ train_and_evaluate_probability_weights(
        partition = prob_partitions[[.x]],
        seed = .x,
        set_tau = tau,
        set_burn_in = burn_in,
        set_iter = iter,
        set_thin = thin,
        return_coverage = compute_coverage,
        set_level = coverage_level
      ),
      .options = furrr_options(seed = TRUE)
    ) |>
      transpose()

    csv_file_name <- paste0("low_data_results_dataset", j, "_n", n, ".csv")

    # write mse
    bind_cols(
      baseline_results$mse |> bind_rows(),
      llm_methods_results$mse |> bind_rows() |> select(-y),
      naive_results$mse |> bind_rows() |> select(-y)
    ) |>
      mutate(
        dataset = j,
        n_requested = n,
        n_realized = realized_n,
        tau = tau,
        .before = 1
      ) |>
      write_csv(csv_file_name)

    # write coverage (only when the samplers retained their draws)
    if (compute_coverage) {
      ss_coverage <- bind_rows(
        baseline_results$coverage |> bind_rows(.id = "partition"),
        llm_methods_results$coverage |> bind_rows(.id = "partition"),
        naive_results$coverage |> bind_rows(.id = "partition")
      )

      if (nrow(ss_coverage) > 0) {
        oos_coverage(ss_coverage, level = coverage_level) |>
          mutate(dataset = j, n_requested = n, tau = tau, .before = 1) |>
          write_csv(paste0("low_data_coverage_dataset", j, "_n", n, ".csv"))
      } else {
        message("    no coverage rows returned; skipping coverage file")
      }
    }

    message("  Completed dataset ", j, ", n = ", n)
  }

  message("Completed dataset ", j)
}

message("low data analysis complete")
