library(rstan)
library(dplyr)

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
  iter = 50000,
  # warmup = 1000,
  chain = 5
)

prior_samples <- extract(model_fit_prior, pars = pars_of_interest, permuted = FALSE) 
index_vec <- apply(prior_samples, 1, function(x_mat) {
  all(!is.na(x_mat))
})
prior_samples <- prior_samples[index_vec, , ,drop = FALSE]

n_prior_samples <- 5000

prior_samples <- prior_samples[
  sample(x = 1 : sum(index_vec), n_prior_samples / 5, replace = T),
  , ,
  drop = FALSE
]

saveRDS(
  object = prior_samples,
  file = "rds/hiv-example/prior-samples.rds"
)
