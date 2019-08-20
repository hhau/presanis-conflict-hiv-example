library(rstan)

data <- readRDS("rds/hiv-example/hiv-data.rds")
prefit <- stan_model("scripts/stan-files/hiv-ev-sythn-small-submodel.stan")
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
  cores = 4,
  chains = 4
)

mcmc_samples <- extract(model_fit, pars = "p_12", permuted = FALSE)
temp_names <- dimnames(mcmc_samples)
temp_names$parameters <- "p[12]"
dimnames(mcmc_samples) <- temp_names

saveRDS(
  object = mcmc_samples,
  file = "rds/hiv-example/small-submodel-samples.rds"
)
