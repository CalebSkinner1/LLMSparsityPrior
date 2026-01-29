# LLM Sparsity Prior for Robust Feature Selection

This repository contains implementation of the **LLM Sparsity Prior**, a Bayesian variable selection framework that integrates *a priori* feature importance synthesized by Large Language Models (LLM) into a Spike-and-Slab prior.

This repository is organized to support three primary objectives:
1) Core implementation: An efficient MCMC Sampler for the LSP
2) Simulations: Code to replicate the comparative simulation studies presented in the paper.
3) Data Application: Prompt engineering templates and analysis code for the Acute Kidney Injury (AKI) application.

## Repository Structure
```{text}
├── MCMC Samplers/
│   ├── LSP_regression_fixed_s.R             # MCMC Sampler for Fixed Sparsity
│   └── LSP_regression_random_s.R            # MCMC Sampler for Random Sparsity
├── Simulations/
│   ├── weight_quality_sims.R     # Main script to run simulations
│   └── weight_quality_support.R  # Support functions for weight generation, metric computation, etc.
├── AKI Data Application/
│   ├── Prompt_Engineering.ipynb    # Exact prompts used for GPT-4o
│   └── End_to_End_Analysis.ipynb   # Preprocessing and model fitting pipeline
└── README.md
```

## Installation and Dependencies
The MCMC samplers and simulations are written in R. To run these, please be sure the following packages are installed:
```{r}
install.packages(c("MASS", "tidyverse", "hypergeo", "glmnet", "parallel", "furrr"))
```

## Usage Example
```{r}
source("LSP_regression_fixed_s.R")

# Generate synthetic data
n <- 100; p <- 100; signals <- 5

X <- MASS::mvrnorm(n, mu = rep(0, p), diag(p))
beta_true <- c(rep(0, p - s), rep(1, s))
alpha_true <- 1
y <- X %*% beta_true + alpha_true + rnorm(n, 0, sd = 1)

# Define LLM-generated feature weights (in practice, these come from the LLM prompt)
weights <- sample(1:5, p, replace = TRUE)

# Run Sampler
fit <- lsp_fixed_gibbs_sampler(
  X = X,
  y = y,
  weights = weights,
  c = 1,
  sparsity = 0.05,
  a_sigma = 1,
  a_sigma = 1,
  tau = 1)

# Extract Posterior Means, etc.
beta_est <- colMeans(fit$beta)
```


