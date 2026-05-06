# LLM Sparsity Prior for Robust Feature Selection

This repository contains an implementation of the **LLM Sparsity Prior**, a Bayesian variable selection framework that integrates *a priori* feature importance synthesized by Large Language Models (LLMs) into the Spike-and-Slab and Spike-and-Slab Lasso priors.

This repository is organized to support three primary objectives:
1) Core implementation: Efficient posterior estimation for LSP methods.
2) Simulations: Comparative simulation studies.
3) Data Application: Prompt engineering templates and Acute Kidney Injury (AKI) analysis.

## Repository Structure
```text
.
├── LSP_SS/
│   ├── LSP_SSR_fixed_s.R             # MCMC Sampler for LSP (SS) with Fixed Sparsity
│   └── LSP_SSR_random_s.R            # MCMC Sampler for LSP (SS) with Random Sparsity
│── LSP_SSL/
│   ├── LSP_SSLR.R                    # R Function for MAP Estimation of LSP (SSL)
│   ├── LSP_SSL_descent.c             # C code for coordinate descent algorithm (lightly edited from https://github.com/cran/SSLASSO)
│   ├── LSP_SSL_functions.c           # C functions for the descent algorithm
├── Simulations/
│   ├── weight_quality_sims.R         # Main script to run simulations
│   ├── eta_sensitivity_sims.R        # Script to run simulations on eta sensitivity analysis
│   └── weight_quality_support.R      # Support functions for weight generation, metric computation, etc.
├── AKI Data Application/
│   ├── analysis/
│   ├──── aki_analysis_support.R      # Support functions for AKI Analysis
│   ├──── aki_subsets.R               # Script to run analysis over five subsets
│   ├──── aki_low_data.R              # Script to run analysis over low-data regime
│   ├──── aki_weight_sensitivity.R    # Script to run sensitivity analysis
│   ├── weight_prompts/
│   ├──── prompty_adjust_collinearity_constraints.ipynb   # Prompt to generate LLM weights from GPT 5.2o (Collinearity Constraints Adjusted)
│   ├──── prompty_adjust_reasoning_structure.ipynb        # Prompt to generate LLM weights from GPT 5.2o (Reasoning Structure Adjusted)
│   ├──── prompty_adjust_scoring_rubric10.ipynb           # Prompt to generate LLM weights from GPT 5.2o (Scoring Rubric Adjusted)
│   ├──── prompty_adjust_task.ipynb                       # Prompt to generate LLM weights from GPT 5.2o (Task Definition Adjusted)
│   ├──── prompty_original.ipynb                          # Prompt to generate LLM weights from GPT 5.2o (Original)
│   ├──── prompty_adjust_probabilities.ipynb              # Prompt to generate Naive LLM weights from GPT 5.2o (Probabilities)
│   ├── weights/                       # directory of generated weights: five prompt variants by five runs and naive weights
└── README.md
```

## Installation and Dependencies
The core method, simulation scripts, and data analysis are written in R. To run these, be sure the following packages are installed:
```r
install.packages(c("MASS", "tidyverse", "hypergeo", "glmnet", "furrr", "pROC", "tidymodels", "simstudy", "Mhorseshoe"))
```

Spike-and-Slab Lasso files must be compiled locally. Run the following command in your terminal at the repository root.
```bash
cd LSP_SSL
R CMD SHLIB LSP_SSL_descent.c LSP_SSL_functions.c -o lsp_ssl.so
```

## Usage Example
```r
source("LSP_SS/LSP_SSR_random_s.R")
source("LSP_SSL/LSP_SSLR.R")

# Generate synthetic data
set.seed(1); n <- 50; p <- 100; signals <- 5

X <- MASS::mvrnorm(n, mu = rep(0, p), diag(p))
beta_true <- c(rep(0, p - signals), rep(1, signals))
alpha_true <- 1
y <- X %*% beta_true + alpha_true + rnorm(n, 0, sd = 1)

# Define LLM-generated feature weights (in practice, these come from the LLM prompt)
weights <- c(rep(1, p - signals), rep(5, signals)) # perfect weights

# Run Sampler with default settings
lsp_ss_fit <- lsp_random_ss_gibbs_sampler(
  X = X,
  y = y,
  weights = weights)

# Extract Posterior Means, etc.
lsp_ss_beta_est <- colMeans(lsp_ss_fit$beta)

lsp_ssl_fit <- lsp_ssl_map(
  X = X,
  y = y,
  weights = weights,
  penalty = "adaptive")

 # Select model from descent (use BIC to select lambda_0)
 lsp_ssl_beta_est <- select_lambda0_bic(lsp_ssl_fit, X = X, y = y)

 # print estimates
 round(lsp_ss_beta_est, digits = 2)
 # [1] 0.76 0.00 0.00 ... 0.00 0.85 1.05 1.11 0.86 1.15
 round(as.vector(lsp_ssl_beta_est), digits = 2)
 # [1] 0.78 0.00 0.00 ... 0.00 0.84 1.05 1.12 0.86 1.15
```