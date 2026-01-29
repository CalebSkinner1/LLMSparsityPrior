# LLM Sparsity Prior for Robust Feature Selection

This repository contains implementation of the **LLM Sparsity Prior**, a Bayesian variable selection framework that integrates *a priori* feature importance synthesized by Large Language Models (LLM) into a Spike-and-Slab prior.

This repository is organized to support three primary objectives:
1) Core implementation: An efficient MCMC Sampler for the LSP
2) Simulations: Code to replicate the comparative simulation studies presented in the paper.
3) Data Application: Prompt engineering templates and analysis code for the Acute Kidney Injury (AKI) application.

## Installation and Dependencies
The MCMC samplers are written in R and require no packages beyond base R. To run the simulations, please be sure the following packages are installed:
```{r}
install.packages(c("tidyverse", "hypergeo", "glmnet", "parallel", "furrr")
```


## MCMC Samplers
We provide code for both the fixed and random sparsity variations of the LLM Sparsity Prior. The MCMC Sampler for fixed sparsity is a Gibbs Sampler with a single Add-Delete-Swap (Metropolis Hastings) step for model selection.
The random sparsity is similar, but with an added Metropolis Hastings Step for estimating the sparsity parameter $s$.

## Simulations
Simulations may be generated at weight_quality_sims.R with parameters adjusted at the top of the file. Helper functions are located in weight_quality_support.R to faciliate readability.

## AKI
The AKI data application necessitates two important files: Prompt_Engineering and End_to_End_Analysis. Prompt Engineering details the specific phrases used to prompt the LLM,
while the End to End Analysis contains data preprocessing and analysis for the paper.


