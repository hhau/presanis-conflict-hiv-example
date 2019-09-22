library(rstan)
library(wsre)
library(ggplot2)

pre_wsre_fit <- stan_model("scripts/stan-files/hiv-ev-sythn-prior-wsre.stan")

wsre_est <- wsre(
  stanmodel = pre_wsre_fit,
  wf_mean = array(c(wf_mean = seq(from = 0.01, to = 0.99, length.out = 10))),
  wf_pars = list(wf_sd = array(c(wf_sd = 0.25)), wf_exponent = 1, target_dimension = 1),
  flog_threshold = futile.logger::TRACE,
  stan_control_params = list(adapt_delta = 0.999, max_treedepth = 14),
  n_mcmc_samples = 250
)

# some_values <- seq(from = 0, to = 1, length.out = 100)
# wsre_values <- evaluate(wsre_est, some_values, some_values)
# plot_df <- expand.grid(x = some_values, y = some_values)
# plot_df$val <- wsre_values

# ggplot(plot_df, aes(x = x, y = y, fill = log(val))) +
#   geom_raster()

saveRDS(
  object = wsre_est,
  file = "rds/hiv-example/big-sub-prior-wsre-est.rds"
)
