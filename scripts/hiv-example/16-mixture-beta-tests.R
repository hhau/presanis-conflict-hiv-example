ref_samples <- readRDS(
  "rds/hiv-example/reference-prior-samples.rds"
)

ref_phi_samples <- ref_samples[, , "p[12]"] %>% as.numeric()
hist(ref_phi_samples, breaks = 50, freq = FALSE)
curve_f <- function(x) dbeta(x, 3.4520971, 0.8341708)
curve(curve_f, add = T)

## compare
qqplot(
  x = ref_phi_samples,
  y = rbeta(n = length(ref_phi_samples * 10), 3.4520971, 0.8341708)
)
abline(a = 0, b = 1)

library(BetaMixture)
res <- BM_Fit(ref_phi_samples, K = 2)

new_curve_f <- function(x) {
  0.4999985 * dbeta(x, 3.299886, 0.8451038) + 
  0.5000015 * dbeta(x, 3.766859, 0.8457738)
}

curve(new_curve_f, add = T, col = "red")

beta_mean <- function(a, b) {
  a / (a + b)
}

beta_variance <- function(a, b) {
  (a * b) / ((a + b)^2 * (a + b + 1))
}

beta_mean(3.299886, 0.8451038)
beta_mean(3.766859, 0.8457738)

beta_variance(3.299886, 0.8451038)
beta_variance(3.766859, 0.8457738)


library(rstan)

stan_prefit <- stan_model(
  model_code = 
  "data {
  } 
  parameters { 
    real <lower = 0, upper = 1> x;
  } 
  model {
    target += log_mix(
      0.4999985, 
      beta_lpdf(x | 3.299886, 0.8451038), 
      beta_lpdf(x | 3.766859, 0.8457738)
    );
  }"
)

stan_samples <- sampling(
  stan_prefit,
  chains = 6,
  iter = 1e5 + 500,
  warmup = 500,
  cores = 6
)

mix_samples <- extract(stan_samples, permuted = FALSE, pars = "x") %>% as.numeric()

qqplot(
  ref_phi_samples,
  mix_samples
)
abline(a = 0, b = 1, col = "red")
