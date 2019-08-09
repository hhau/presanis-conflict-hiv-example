source("scripts/common/plot-settings.R")

library(bayesplot)
library(ggplot2)

prior_samples <- readRDS(file = "rds/hiv-example/prior-samples.rds")
full_model_sample <- readRDS(file = "rds/hiv-example/full-model-fit.rds")
big_submodel_sample <- readRDS(file = "rds/hiv-example/big-submodel-samples.rds")
pars_of_interest <- readRDS(file = "rds/hiv-example/pars-of-interest.rds")

# Prior / posterior mcmc-areas in red/blue to show concentration?
bp_prior <- bayesplot::mcmc_intervals(
  prior_samples,
  pars = pars_of_interest
)

bp_subpost <- bayesplot::mcmc_intervals(
  big_submodel_sample,
  pars = pars_of_interest
)

bp_post <- bayesplot::mcmc_intervals(
  full_model_sample,
  pars = pars_of_interest
)

prior_data <- bp_prior$data
subpost_data <- bp_subpost$data
post_data <- bp_post$data

prior_data$dtype <- "cc_prior"
subpost_data$dtype <- "bb_subpost"
post_data$dtype <- "aa_post"

full_data <- dplyr::bind_rows(prior_data, subpost_data, post_data)

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
      cc_prior = as.character(blues[2]),
      bb_subpost = as.character(greens[2]),
      aa_post = highlight_col
    ),
    labels = c(
      cc_prior = expression("p"(phi)),
      bb_subpost = expression("p"[2](phi~"|"~"Y"[2])),
      aa_post = expression("p"(phi~"|"~"Y"[1]~","~"Y"[2]))
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
