source("scripts/common/plot-settings.R")

library(bayesplot)
library(ggplot2)
library(tidyr)
library(dplyr)
library(abind)

prior_samples <- readRDS(file = "rds/hiv-example/prior-samples.rds")
full_model_sample <- readRDS(file = "rds/hiv-example/full-model-fit.rds")
big_submodel_sample <- readRDS(file = "rds/hiv-example/big-submodel-samples.rds")
small_submodel_sample <- readRDS(file = "rds/hiv-example/small-submodel-samples.rds")

# stage one here is targeting the big submodel with the 
# prior on phi marginalised out
stage_one_samples <- readRDS(file = "rds/hiv-example/stage-one-samples.rds")

# wsre'd stage one samples
wsre_stage_one_samples <- readRDS(file = "rds/hiv-example/stage-one-wsre-samples.rds")
## TODO:
# - add p12, compute it and dimnames it - abind it onto the rest of the samples
res <- apply(wsre_stage_one_samples, 1:2, function(x) {
  db <- x["dd"] * x["bb"]
  e1ab <- x["ee"] * (1 - x["aa"] - x["bb"])
  p12_prop <<- (db + x["ww"] * e1ab) / (db + e1ab)
  return(p12_prop)
})
n_iter <- dim(wsre_stage_one_samples)[1]
n_chain <- dim(wsre_stage_one_samples)[2]
res <- array(as.vector(res), dim = c(n_iter, n_chain, 1))
old_dimnames <- dimnames(wsre_stage_one_samples)
old_dimnames[[3]] <- c(old_dimnames[[3]], "p[12]")
wsre_stage_one_samples <- abind(wsre_stage_one_samples, res, along = 3)
dimnames(wsre_stage_one_samples) <- old_dimnames

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

bp_small_subpost <- bayesplot::mcmc_intervals(
  small_submodel_sample
)

bp_post <- bayesplot::mcmc_intervals(
  full_model_sample,
  pars = pars_of_interest
)

bp_stage_one <- bayesplot::mcmc_intervals(
  stage_one_samples,
  pars = pars_of_interest
)

bp_wsre_stage_one <- bayesplot::mcmc_intervals(
  wsre_stage_one_samples
)

prior_data <- bp_prior$data
subpost_data <- bp_subpost$data
small_subpost_data <- bp_small_subpost$data
post_data <- bp_post$data
stage_one_data <- bp_stage_one$data
wsre_stage_one_data <- bp_wsre_stage_one$data

prior_data$dtype <- "cc_prior"
subpost_data$dtype <- "bb_subpost"
small_subpost_data$dtype = "ba_subpost"
wsre_stage_one_data$dtype <- "ab_wsre_stage_one_target"
stage_one_data$dtype <- "ab_stage_one_target"
post_data$dtype <- "aa_post"

full_data <- dplyr::bind_rows(prior_data, subpost_data, small_subpost_data, post_data, stage_one_data, wsre_stage_one_data)
full_data <- tidyr::complete(full_data, dtype, parameter)

## TODO: move this else where, need more colours for the moment
purples <- RColorBrewer::brewer.pal(4, name = "Purples")

# We need to do some tidy-ing here
# The p[\d+]'s need zero-padding.
# The base parameters need prepending with base_
# zero padding first 
values <- sprintf("p[%02d]", 1:12)
names(values) <- sprintf("p[%d]", 1:12)
full_data$parameter <- full_data$parameter %>% 
  recode(!!!values)

# now the base prepending
prepend_values <- sapply(
  X = c("aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh", "ww"),
  function(x) {
    paste("base_", x, sep = "")
  }
) 
full_data$parameter <- full_data$parameter %>% 
  recode(!!!prepend_values)

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
      bb_subpost = as.character(greens[1]),
      ba_subpost = as.character(greens[4]),
      ab_wsre_stage_one_target = purples[2],
      ab_stage_one_target = purples[4],
      aa_post = highlight_col
    ),
    labels = c(
      cc_prior = expression("p"(phi)),
      bb_subpost = expression("p"[2](phi~"|"~"Y"[2])),
      ba_subpost = expression("p"[1](phi~"|"~"Y"[1])),
      ab_wsre_stage_one_target = expression(hat("p")[2](phi~"|"~"Y"[2])~"/"~"p"[2](phi)),
      ab_stage_one_target = expression("p"[2](phi~"|"~"Y"[2])~"/"~"p"[2](phi)),
      aa_post = expression("p"(phi~"|"~"Y"))
    ),
    guide = guide_legend(reverse = TRUE)
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
