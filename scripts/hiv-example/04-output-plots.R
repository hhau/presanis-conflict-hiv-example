source("scripts/common/plot-settings.R")

library(bayesplot)
library(ggplot2)

prior_samples <- readRDS(file = "rds/hiv-example/prior-samples.rds")
full_model_sample <- readRDS(file = "rds/hiv-example/full-model-fit.rds")
pars_of_interest <- readRDS(file = "rds/hiv-example/pars-of-interest.rds")

# Prior / posterior mcmc-areas in red/blue to show concentration?
bp_prior <- bayesplot::mcmc_intervals(
  prior_samples,
  pars = pars_of_interest
)

bp_post <- bayesplot::mcmc_intervals(
  full_model_sample,
  pars = pars_of_interest
)

prior_data <- bp_prior$data
post_data <- bp_post$data

prior_data$dtype <- "prior"
post_data$dtype <- "post"

full_data <- dplyr::bind_rows(prior_data, post_data)

p1 <- ggplot(full_data, aes(x = parameter, group = interaction(parameter, dtype), col = dtype)) +
  geom_boxplot(
    aes(
      ymin = ll,
      lower = l,
      middle = m,
      upper = h,
      ymax = hh
    ),
    stat = "identity"
  ) +
  scale_discrete_manual(
    aesthetics = "col",
    values = c(
      prior = as.character(blues[2]),
      post = highlight_col
    ),
    labels = c(
      prior = "Prior",
      post = "Posterior"
    )
  ) +
  labs(col = "QoI") +
  xlab("Parameter") +
  ylab("Probability") +
  coord_flip() +
  NULL

ggsave_fullpage(
  filename = "plots/hiv-example/prior-post-compare.pdf",
  plot = p1
)
