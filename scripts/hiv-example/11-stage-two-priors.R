# This file is made to be source'd into the stage two melded posterior file.
library(wsre)
library(futile.logger)

# set up the functions
# big submodel wrse - needed if we use logarithmic pooling
flog.info("Reading samples and wsre estimate")
big_model_wsre <- readRDS(
  file = "rds/hiv-example/big-sub-prior-wsre-est.rds"
)

# big submodel prior  samples
big_model_prior_samples <- readRDS(file = "rds/hiv-example/prior-samples.rds")
p12_big_model_samples <- as.matrix(as.vector(big_model_prior_samples[, 1:5, "p[12]"]))
p12_bw <- bw.SJ(p12_big_model_samples)

flog.info("Building stage two prior marginals")
log_big_p12_prior_marginal <- function(x) {
  log(wsre:::kde_func_nd(x, p12_big_model_samples, p12_bw))
}

log_big_p12_prior_marginal_beta_approx <- function(x) {
  # magic numbers from MLE using 100,000 samples from the prior
  # visually gives very good fit.
  dbeta(x, 3.4520971, 0.8341708, log = TRUE)
}

# stage two / Small submodel prior marginal (free to choose this)
log_small_p12_prior_marginal <- function(x) {
  dbeta(x, 1, 1, log = TRUE)
}
  
# pooled prior
log_logarithmic_pooled_prior_no_wsre <- function(phi_nu, phi_de) {
  0.5 * (
    log_big_p12_prior_marginal(phi_nu) -
    log_big_p12_prior_marginal(phi_de)
  )
}

log_logarithmic_pooled_prior_with_wsre <- function(phi_nu, phi_de) {
  0.5 * log(evaluate(big_model_wsre, phi_nu, phi_de))
}

log_logarithmic_pooled_prior_with_wsre_tele <- function(phi_nu, phi_de) {
  0.5 * log(evaluate_telescope_fixed_dist(big_model_wsre, phi_nu, phi_de))
}

log_logarithmic_pooled_prior_beta_approx <- function(phi_nu, phi_de) {
  0.5 * (
    log_big_p12_prior_marginal_beta_approx(phi_nu) -
    log_big_p12_prior_marginal_beta_approx(phi_de)
  )
}

# should really LSE this
# there is no way to use the wsre estimate in the pooled prior
# for the linear pooling case
log_linear_pooled_prior_no_wsre <- function(phi_nu, phi_de) {
  log(
    exp(log_p12_prior_marginal(phi_nu))
    +
    exp(log_small_p12_prior_marginal(phi_nu))
  )
  -
  log(
    exp(log_p12_prior_marginal(phi_de))
    +
    exp(log_small_p12_prior_marginal(phi_de))
  )
}
