library(tibble)
library(dplyr)
library(wsre)
library(pbapply)

source("scripts/hiv-example/18-telescoping-tests.R")

ref_stage_one_samples <- readRDS(
  "rds/hiv-example/stage-one-reference-samples.rds"
)

ref_melded_posterior_samples <- readRDS(
  "rds/hiv-example/stage-two-reference-samples.rds"
)

wsre_estimate <- readRDS(
  "rds/hiv-example/big-sub-prior-wsre-est.rds"
)

naive_prior_samples <- readRDS(
  file = "rds/hiv-example/prior-samples.rds"
)

n_estimates <- length(wsre_estimate$estimates)
no_naive_wsre_estimate <- wsre_estimate
no_naive_wsre_estimate$estimates[[n_estimates]] <- NULL

naive_kde_bw <- bw.SJ(naive_prior_samples[, , "p[12]"] %>% as.numeric())
naive_kde <- function(x) {
  wsre:::kde_func_nd(
    x, 
    naive_prior_samples[, , "p[12]"] %>% 
      as.numeric() %>% 
      as.matrix()
    , naive_kde_bw
  )
}

naive_ratio <- function(x_nu, x_de) {
  naive_kde(x_nu) / naive_kde(x_de)
}

ref_ratio <- function(x_nu, x_de) {
  dbeta(x_nu, 3.4520971, 0.8341708) / dbeta(x_de, 3.4520971, 0.8341708) 
}

# set up some denominator values
x_de_vals <- seq(
  from = 0.15, 
  to = 0.75, 
  length.out = 9
)

# setup up a fixed set 
x_nu_vals <- seq(from = 0, to = 1, length.out = 250)

# takes about ~15 mins
res <- pblapply(x_de_vals, cl = 6, function(a_de_val) {
  tibble(
    x_nu = rep(x_nu_vals, 9),
    x_de = rep(a_de_val, 9 * length(x_nu_vals)),
    value = c(
      ref_ratio(x_nu_vals, a_de_val),
      sapply(x_nu_vals, function(a_nu_val) naive_ratio(a_nu_val, a_de_val)),
      sapply(x_nu_vals, function(a_nu_val) evaluate(wsre_estimate, a_nu_val, a_de_val, mc_cores = 1)),
      sapply(x_nu_vals, function(a_nu_val) evaluate_telescope_1d(wsre_estimate, a_nu_val, a_de_val)),
      sapply(x_nu_vals, function(a_nu_val) evaluate_telescope_1d(no_naive_wsre_estimate, a_nu_val, a_de_val)),
      sapply(x_nu_vals, function(a_nu_val) evaluate_telescope_1d_means(no_naive_wsre_estimate, a_nu_val, a_de_val)),
      sapply(x_nu_vals, function(a_nu_val) evaluate_telescope_1d_dumb(wsre_estimate, a_nu_val, a_de_val)),
      sapply(x_nu_vals, function(a_nu_val) evaluate_telescope_fixed_N(wsre_estimate, a_nu_val, a_de_val)),
      sapply(x_nu_vals, function(a_nu_val) evaluate_telescope_fixed_dist(wsre_estimate, a_nu_val, a_de_val))
    ),
    d_type = rep(
      c(
        "ref", 
        "naive", 
        "wsre", 
        "wsre_tele1d", 
        "wsre_tele1d_no_naive", 
        "wsre_tele1d_no_naive_means", 
        "wsre_tele_dumb",
        "wsre_tele_fixed_N",
        "wsre_tele_fixed_dist"
      ), 
      each = length(x_nu_vals)
    )
  )
}) 

plot_df <- do.call(rbind, res)
plot_df$x_de_plot <- sprintf("'x'['de'] ~ '=' ~ %.2f", plot_df$x_de)

saveRDS(
  file = "rds/hiv-example/ratio-estimates-df.rds",
  object = plot_df
)

# do the inverse
x_nu_vals_inverse <- seq(
  from = 0.15, 
  to = 0.75, 
  length.out = 9
)

# setup up a fixed set 
x_de_vals_inverse <- seq(from = 0, to = 1, length.out = 250)

# takes about ~15 mins
res_inverse <- pblapply(x_nu_vals_inverse, cl = 6, function(a_nu_val) {
  tibble(
    x_de = rep(x_de_vals_inverse, 9),
    x_nu = rep(a_nu_val, 9 * length(x_de_vals_inverse)),
    value = c(
      ref_ratio(a_nu_val, x_de_vals_inverse),
      sapply(x_de_vals_inverse, function(a_de_val) naive_ratio(a_nu_val, a_de_val)),
      sapply(x_de_vals_inverse, function(a_de_val) evaluate(wsre_estimate, a_nu_val, a_de_val, mc_cores = 1)),
      sapply(x_de_vals_inverse, function(a_de_val) evaluate_telescope_1d(wsre_estimate, a_nu_val, a_de_val)),
      sapply(x_de_vals_inverse, function(a_de_val) evaluate_telescope_1d(no_naive_wsre_estimate, a_nu_val, a_de_val)),
      sapply(x_de_vals_inverse, function(a_de_val) evaluate_telescope_1d_means(no_naive_wsre_estimate, a_nu_val, a_de_val)),
      sapply(x_de_vals_inverse, function(a_de_val) evaluate_telescope_1d_dumb(wsre_estimate, a_nu_val, a_de_val)),
      sapply(x_de_vals_inverse, function(a_de_val) evaluate_telescope_fixed_N(wsre_estimate, a_nu_val, a_de_val)),
      sapply(x_de_vals_inverse, function(a_de_val) evaluate_telescope_fixed_dist(wsre_estimate, a_nu_val, a_de_val))
    ),
    d_type = rep(
      c(
        "ref", 
        "naive", 
        "wsre", 
        "wsre_tele1d", 
        "wsre_tele1d_no_naive", 
        "wsre_tele1d_no_naive_means", 
        "wsre_tele_dumb",
        "wsre_tele_fixed_N",
        "wsre_tele_fixed_dist"
      ), 
      each = length(x_de_vals_inverse)
    )
  )
}) 

plot_df_inverse <- do.call(rbind, res_inverse)
plot_df_inverse$x_nu_plot <- sprintf("'x'['nu'] ~ '=' ~ %.2f", plot_df_inverse$x_nu)

saveRDS(
  file = "rds/hiv-example/ratio-estimates-df-inverse.rds",
  object = plot_df_inverse
)
