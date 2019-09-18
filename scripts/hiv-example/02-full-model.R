library(rstan)

data <- readRDS("rds/hiv-example/hiv-data.rds")
prefit <- stan_model("scripts/stan-files/hiv-ev-sythn.stan")
pars_of_interest <- readRDS(
  file = "rds/hiv-example/pars-of-interest.rds"
)

stan_data <- list(
  n_studies = nrow(data),
  y_obs = data$y,
  n_obs = data$n
)

model_fit <- sampling(
  prefit,
  data = stan_data,
  cores = 5,
  chains = 5,
  iter = 2000,
  warmup = 1000,
  control = list(adapt_delta = 0.9, max_treedepth = 12)
)

mcmc_samples <- extract(model_fit, pars = pars_of_interest, permuted = FALSE)

saveRDS(
  object = mcmc_samples,
  file = "rds/hiv-example/full-model-fit.rds"
)


