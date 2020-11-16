library(tibble)
library(dplyr)
library(wsre)
library(pbapply)

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
  length.out = 6
)

# setup up a fixed set 
x_nu_vals <- seq(from = 0, to = 1, length.out = 250)

# takes about ~15 mins
res <- pblapply(x_de_vals, cl = 6, function(a_de_val) {
  tibble(
    x_nu = rep(x_nu_vals, 3),
    x_de = rep(a_de_val, 3 * length(x_nu_vals)),
    value = c(
      ref_ratio(x_nu_vals, a_de_val),
      sapply(x_nu_vals, function(a_nu_val) naive_ratio(a_nu_val, a_de_val)),
      sapply(x_nu_vals, function(a_nu_val) evaluate(wsre_estimate, a_nu_val, a_de_val))
    ),
    d_type = rep(c("ref", "naive", "wsre"), each = length(x_nu_vals))
  )
}) 

plot_df <- do.call(rbind, res)
plot_df$x_de_plot <- sprintf("'x'['de'] ~ '=' ~ %.2f", plot_df$x_de)

saveRDS(
  file = "rds/hiv-example/ratio-estimates-df.rds",
  object = plot_df
)
