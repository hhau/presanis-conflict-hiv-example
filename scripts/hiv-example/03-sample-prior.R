library(rstan)

data <- readRDS("rds/hiv-example/hiv-data.rds")
prefit <- stan_model("scripts/stan-files/hiv-ev-sythn-prior.stan")
pars_of_interest <- readRDS(
  file = "rds/hiv-example/pars-of-interest.rds"
)

## prior samples
stan_data_prior <- list(
  n_studies = nrow(data)
)

model_fit_prior <- sampling(
  prefit,
  data = stan_data_prior,
  algorithm = "Fixed_param",
  iter = 20000,
  chain = 5
)

prior_samples <- extract(model_fit_prior, pars = pars_of_interest, permuted = FALSE)

saveRDS(
  object = prior_samples,
  file = "rds/hiv-example/prior-samples.rds"
)
