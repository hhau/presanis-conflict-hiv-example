library(rstan)
library(dplyr)

data <- readRDS("rds/hiv-example/hiv-data.rds")
pars_of_interest <- readRDS("rds/hiv-example/pars-of-interest.rds")

prefit <- stan_model(
  "scripts/stan-files/hiv-ev-sythn-stage-one-target-big-reference.stan"
)

stan_data <- list(
  n_studies = nrow(data),
  y_obs = data$y,
  n_obs = data$n
)

model_fit <- sampling(
  prefit,
  data = stan_data,
  cores = 6,
  chains = 6,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.9, max_treedepth = 12)
)

mcmc_samples <- extract(model_fit, pars = pars_of_interest, permuted = FALSE)

saveRDS(
  object = mcmc_samples,
  file = "rds/hiv-example/stage-one-reference-samples.rds"
)