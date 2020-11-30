library(rstan)
library(dplyr)

prefit <- stan_model("scripts/stan-files/hiv-ev-sythn-stage-one-target-big.stan")
prior_samples <- readRDS("rds/hiv-example/prior-samples.rds")

prior_index <- grep("12", dimnames(prior_samples)$parameters)
p12_prior_samples <- prior_samples[, , prior_index] %>% as.vector()

data <- readRDS("rds/hiv-example/hiv-data.rds")
pars_of_interest <- readRDS(
  file = "rds/hiv-example/pars-of-interest.rds"
)

stan_data <- list(
  n_studies = nrow(data),
  y_obs = data$y,
  n_obs = data$n,
  n_prior_samples = length(p12_prior_samples),
  phi_prior_samples = p12_prior_samples,
  bandwidth = bw.SJ(p12_prior_samples)
)

model_fit <- sampling(
  prefit,
  data = stan_data,
  cores = 6,
  chains = 24,
  iter = 20000,
  warmup = 2000,
  control = list(adapt_delta = 0.9, max_treedepth = 12)
)

mcmc_samples <- extract(model_fit, pars = pars_of_interest, permuted = FALSE)

saveRDS(
  object = mcmc_samples,
  file = "rds/hiv-example/stage-one-samples.rds"
)