# AKI Data Application - Single Run
#
# Applies all LSP to five clinical subgroups derived from the AKI dataset.
# Does not perform repeated cross validation; trains on entire dataset

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
# Two weight sets are used:
#   aki_weights_0             — discretized LLM importance weights (for LSP-SS/SSL)
#   aki_weights_probabilities — probability importance weights (used as direct
#                            prior inclusion probabilities)
#
# subset_weights aligns the weight data frame to the columns present in each
# subgroup dataset after topK_features filtering.
# ------------------------------------------------------------------------------

aki_weights_0 <- read_csv(
  "AKI_Data_Application/weights/aki_weights_original_1.csv",
  show_col_types = FALSE
)
aki_weights_probabilities <- read_csv(
  "AKI_Data_Application/weights/aki_weights_probabilities_1.csv",
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

tau <- 2
eta_range <- NULL # NULL triggers default prior
iter <- 125000
burn_in <- 25000
random_s <- TRUE
fixed_s <- FALSE

n_top_models <- 5
credible_level <- 0.95

# ------------------------------------------------------------------------------
# Estimate phi hat
# ------------------------------------------------------------------------------

message("estimating phi hat...")

phi_estimates <- map2(
  list(aki_data1, aki_data2, aki_data3, aki_data4, aki_data5),
  list(aki_weights1, aki_weights2, aki_weights3, aki_weights4, aki_weights5),
  ~ estimate_phi(
    data = .x,
    outcome_var = outcome,
    weights = .y,
    model = "SSL",
    set_tau = tau,
    iter = iter,
    burn_in = burn_in
  )
)

phi_estimates

# ------------------------------------------------------------------------------
# Run full models
# ------------------------------------------------------------------------------

focus_data <- aki_data1
focus_weights <- aki_weights1
message("fitting full models on subgroup: 1")
# may adjust to aki_data2, etc.

scaled_data <- train_test_split(
  focus_data,
  focus_weights,
  outcome_var = outcome,
  n_folds = 1
) |>
  pluck(1)

feature_names <- colnames(scaled_data$X_train_scaled)
importance <- scaled_data$weights$importance

lsp_full <- lsp_random_ss_gibbs_sampler(
  X = scaled_data$X_train_scaled,
  y = scaled_data$y_train_scaled,
  weights = importance,
  tau = tau,
  iter = iter,
  burn_in = burn_in
)

ss_full <- lsp_random_ss_gibbs_sampler(
  X = scaled_data$X_train_scaled,
  y = scaled_data$y_train_scaled,
  weights = NULL,
  tau = tau,
  iter = iter,
  burn_in = burn_in
)

n_draws <- iter - burn_in

# ------------------------------------------------------------------------------
# Posterior of eta
# ------------------------------------------------------------------------------

eta_summary <- tibble(
  mean = mean(lsp_run$eta),
  median = median(lsp_run$eta),
  q10 = quantile(lsp_run$eta, 0.10),
  q90 = quantile(lsp_run$eta, 0.90)
)

eta_summary

eta_plot <- tibble(eta = lsp_run$eta) |>
  ggplot(aes(x = eta)) +
  geom_histogram(bins = 30) +
  labs(x = expression(eta), y = "Posterior draws") +
  theme_minimal()

eta_plot

# Draws where the implied maximum inclusion probability exceeds 1
theta_max_draws <- tibble(eta = lsp_run$eta, s = lsp_run$s) |>
  mutate(
    theta_max = s * map_dbl(eta, \(e) max(importance)^e / mean(importance^e))
  ) |>
  filter(theta_max > 1)

theta_max_draws |> nrow()

# ------------------------------------------------------------------------------
# Highest Posterior models
# ------------------------------------------------------------------------------

# Rank the most frequently visited inclusion patterns in a gamma sample matrix
top_models <- function(gamma, feature_names, n_top = 5) {
  patterns <- apply(gamma, 1, function(g) {
    if (all(g == 0)) {
      "None (all zero)"
    } else {
      paste(feature_names[g == 1], collapse = " + ")
    }
  })

  counts <- sort(table(patterns), decreasing = TRUE)

  tibble(
    combination = names(counts),
    frequency = as.integer(counts)
  ) |>
    slice_head(n = n_top)
}

final_ranking <- top_models(lsp_run$gamma, feature_names, n_top = n_top_models)

final_ranking

# ------------------------------------------------------------------------------
# Marginal inclusion probabilities
# ------------------------------------------------------------------------------

mip <- tibble(
  feature = feature_names,
  lsp_mip = colMeans(lsp_run$gamma),
  ss_mip = colMeans(ss_run$gamma),
  weight = importance
) |>
  arrange(desc(lsp_mip))

mip |>
  arrange(desc(ss_mip))


# predictive intervals

pred_mat <- lsp_run$beta %*%
  t(cbind(1, scaled_data$X_train_scaled)) +
  matrix(
    rnorm(
      n = (iter - burn_in) * length(scaled_data$y_train_scaled),
      mean = 0,
      sd = sqrt(tau / run$invsigma_2)
    ),
    nrow = (iter - burn_in),
    ncol = length(scaled_data$y_train_scaled),
    byrow = FALSE
  )

int <- apply(pred_mat, 2, quantile, probs = c(0.025, 0.975))

between(scaled_data$y_train_scaled, left = int[1, ], right = int[2, ]) |>
  mean()
