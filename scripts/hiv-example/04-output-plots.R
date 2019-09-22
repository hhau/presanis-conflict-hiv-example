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
  p12_prop <- (db + x["ww"] * e1ab) / (db + e1ab)
  return(p12_prop)
})
n_iter <- dim(wsre_stage_one_samples)[1]
n_chain <- dim(wsre_stage_one_samples)[2]
res <- array(as.vector(res), dim = c(n_iter, n_chain, 1))
old_dimnames <- dimnames(wsre_stage_one_samples)
old_dimnames[[3]] <- c(old_dimnames[[3]], "p[12]")
wsre_stage_one_samples <- abind(wsre_stage_one_samples, res, along = 3)
dimnames(wsre_stage_one_samples) <- old_dimnames

# TODO: read in the stage two samples (and the stage two wsre samples!)
stage_two_samples <- readRDS("rds/hiv-example/stage-two-samples.rds")
stage_two_p12_samples <- stage_two_samples$phi_samples
old_s2_dimnames <- dimnames(stage_two_p12_samples)
new_s2_dimnames <- old_s2_dimnames
new_s2_dimnames[[3]] <- "p[12]"
dimnames(stage_two_p12_samples) <- new_s2_dimnames

wsre_stage_two_samples <- readRDS("rds/hiv-example/stage-two-wsre-samples.rds")
wsre_stage_two_p12_samples <- wsre_stage_two_samples$phi_samples
wsre_old_s2_dimnames <- dimnames(wsre_stage_two_p12_samples)
wsre_new_s2_dimnames <- wsre_old_s2_dimnames
wsre_new_s2_dimnames[[3]] <- "p[12]"
dimnames(wsre_stage_two_p12_samples) <- wsre_new_s2_dimnames

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

bp_stage_two <- bayesplot::mcmc_intervals(
  stage_two_p12_samples
)

bp_wsre_stage_two <- bayesplot::mcmc_intervals(
  wsre_stage_two_p12_samples
)

prior_data <- bp_prior$data
subpost_data <- bp_subpost$data
small_subpost_data <- bp_small_subpost$data
post_data <- bp_post$data
stage_one_data <- bp_stage_one$data
wsre_stage_one_data <- bp_wsre_stage_one$data
stage_two_data <- bp_stage_two$data
wsre_stage_two_data <- bp_wsre_stage_two$data


prior_data$dtype <- "h_prior"
subpost_data$dtype <- "g_subpost"
small_subpost_data$dtype = "f_subpost"
wsre_stage_one_data$dtype <- "e_wsre_stage_one_target"
stage_one_data$dtype <- "d_stage_one_target"
wsre_stage_two_data$dtype <- "c_wsre_stage_two_target"
stage_two_data$dtype <- "b_stage_two_target"
post_data$dtype <- "a_post"

full_data <- dplyr::bind_rows(
  prior_data,
  subpost_data,
  small_subpost_data,
  post_data,
  stage_one_data,
  wsre_stage_one_data,
  stage_two_data,
  wsre_stage_two_data
)

full_data <- tidyr::complete(full_data, dtype, parameter)

## TODO: move this else where, need more colours for the moment
purples <- RColorBrewer::brewer.pal(4, name = "Purples")
oranges <- RColorBrewer::brewer.pal(4, name = "Oranges")

# We need to do some tidy-ing here
# The p[\d+]'s need zero-padding.
# The base parameters need prepending with base_
# zero padding first 
values <- sprintf("italic(p)[%02d]", 1:12)
names(values) <- sprintf("p[%d]", 1:12)
full_data$parameter <- full_data$parameter %>% 
  recode(!!!values)

# now the base prepending
prepend_values <- sapply(
  X = c("aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh", "ww"),
  function(x) {
    str <- substr(x, 1, 1)
    res <- paste0('italic(', str, ')')
  }
) 
full_data$parameter <- full_data$parameter %>% 
  recode(!!!prepend_values) %>% 
  as.factor()

# Hack because I want `w` to come before the probabilities
limits_vec <- levels(full_data$parameter)
new_limits_vec <- c(
  limits_vec[1 : 8],
  limits_vec[21],
  limits_vec[9 : 20]
)

plot_pars <- c(
  "italic(p)[12]",
  "italic(a)",
  "italic(b)",
  "italic(d)",
  "italic(e)",
  "italic(w)"
)

sub_data <- full_data %>% 
  filter(parameter %in% plot_pars)


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
      h_prior = as.character(blues[2]),
      g_subpost = as.character(greens[1]),
      f_subpost = as.character(greens[4]),
      e_wsre_stage_one_target = purples[2],
      d_stage_one_target = purples[4],
      c_wsre_stage_two_target = oranges[2],
      b_stage_two_target = oranges[4],
      a_post = highlight_col
    ),
    labels = c(
      h_prior = expression("p"[2](phi)),
      g_subpost = expression("p"[2](phi~"|"~"Y"[2])),
      f_subpost = expression("p"[1](phi~"|"~"Y"[1])),
      e_wsre_stage_one_target = expression("p"[2](phi~"|"~"Y"[2])~"/"~hat("p'")[2](phi)),
      d_stage_one_target = expression("p"[2](phi~"|"~"Y"[2])~"/"~hat("p")[2](phi)),
      c_wsre_stage_two_target = expression(hat("p'")["meld"](phi~"|"~"Y"[1],~"Y"[2])),
      b_stage_two_target = expression(hat("p")["meld"](phi~"|"~"Y"[1],~"Y"[2])),
      a_post = expression("p"(phi~"|"~"Y"))
    ),
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(col = "QoI") +
  xlab("Parameter") +
  ylab("Probability") +
  scale_x_discrete(
    limits = new_limits_vec,
    labels = function(x) parse(text = x)
  ) +
  coord_flip() +
  NULL

p2 <- ggplot(sub_data, aes(x = parameter, group = interaction(parameter, dtype), col = dtype)) +
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
      h_prior = as.character(blues[2]),
      g_subpost = as.character(greens[1]),
      f_subpost = as.character(greens[4]),
      e_wsre_stage_one_target = purples[2],
      d_stage_one_target = purples[4],
      c_wsre_stage_two_target = oranges[2],
      b_stage_two_target = oranges[4],
      a_post = highlight_col
    ),
    labels = c(
      h_prior = expression("p"[2](phi)),
      g_subpost = expression("p"[2](phi~"|"~"Y"[2])),
      f_subpost = expression("p"[1](phi~"|"~"Y"[1])),
      e_wsre_stage_one_target = expression("p"[2](phi~"|"~"Y"[2])~"/"~hat("p'")[2](phi)),
      d_stage_one_target = expression("p"[2](phi~"|"~"Y"[2])~"/"~hat("p")[2](phi)),
      c_wsre_stage_two_target = expression(hat("p'")["meld"](phi~"|"~"Y"[1],~"Y"[2])),
      b_stage_two_target = expression(hat("p")["meld"](phi~"|"~"Y"[1],~"Y"[2])),
      a_post = expression("p"(phi~"|"~"Y"))
    ),
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(col = "QoI") +
  xlab("Parameter") +
  ylab("Probability") +
  scale_x_discrete(
    limits = new_limits_vec[new_limits_vec %in% plot_pars],
    labels = function(x) parse(text = x)
  ) +
  coord_flip() +
  NULL

ggsave_fullpage(
  filename = "plots/hiv-example/prior-post-compare.pdf",
  plot = p1
)

ggsave_halfheight(
  filename = "plots/hiv-example/p12-prior-post-compare.pdf",
  plot = p2
)

# New version of plot

# get the p[12] only data
new_plot_tbl <- sub_data %>% 
  filter(parameter == "italic(p)[12]") %>% 
  mutate(parameter = as.character(parameter))

# get the new data for the p[1](phi) prior
p1_phi_prior_samples <- array(
  data = runif(n = 2500),
  dim = c(500, 5, 1),
  dimnames = list(
    sprintf("Iteration: %d", 1:500),
    sprintf("Chain: %d", 1:5),
    "italic(p)[12]"
  )
)

p1_phi_prior_data <- bayesplot::mcmc_intervals_data(
  p1_phi_prior_samples
)
p1_phi_prior_data$dtype <- "hh_prior"
p1_phi_prior_data <- p1_phi_prior_data %>% 
  mutate(parameter = as.character(parameter))

boxplot_gap <- 1
x_axis_data <- tibble(
  dtype = c(
    new_plot_tbl$dtype,
    p1_phi_prior_data$dtype
  ), 
  x_val = c(
    0, # a_post
    1, # b_stage_two_target
    2, # c_wsre_stage_two_target
    3 + boxplot_gap, # d_stage_one_target
    4 + boxplot_gap, # e_wsre_stage_one_target
    5 + 2 * boxplot_gap, # f_subpost
    6 + 2 * boxplot_gap, # g_subpost
    7 + 3 * boxplot_gap, # h_prior
    8 + 3 * boxplot_gap # hh_prior
  ),
  fill = c(
    highlight_col, # a_post
    oranges[3], # b_stage_two_target
    oranges[1], # c_wsre_stage_two_target
    purples[3], # d_stage_one_target
    purples[1], # e_wsre_stage_one_target
    greens[3], # f_subpost
    greens[2], # g_subpost
    blues[3], # h_prior
    blues[2] # hh_prior
  )
)

p12_only_data <- new_plot_tbl %>% 
  bind_rows(p1_phi_prior_data) %>% 
  left_join(x_axis_data, by = "dtype") %>% 
  mutate(x_val = as.integer(x_val))

p3 <- ggplot(p12_only_data, aes(x = x_val, group = dtype)) +
  geom_boxplot(
    aes(
      ymin = ll,
      lower = l,
      middle = m,
      upper = h,
      ymax = hh
    ),
    stat = "identity",
    fill = p12_only_data$fill
  ) +
  coord_flip() + 
  scale_x_discrete(
    limits = p12_only_data$x_val,
    labels = rev(
      c(
        "11" = expression("p"[1](phi)),
        "10" = expression("p"[2](phi)),
        "8" = expression("p"[2](phi~"|"~"Y"[2])),
        "7" = expression("p"[1](phi~"|"~"Y"[1])),
        "5" = expression("p"[2](phi~"|"~"Y"[2])~"/"~hat("p'")[2](phi)),
        "4" = expression("p"[2](phi~"|"~"Y"[2])~"/"~hat("p")[2](phi)),
        "2" = expression(hat("p'")["meld"](phi~"|"~"Y"[1],~"Y"[2])),
        "1" = expression(hat("p")["meld"](phi~"|"~"Y"[1],~"Y"[2])),
        "0" = expression("p"(phi~"|"~"Y"))
      )
    )
  ) +
  xlab("") +
  ylab(expression(phi)) +
  ggtitle(
    expression(phi=="p"[12]~": Melding distributions")
  ) +
  theme(
    axis.text = element_text(size = rel(1.0))
  ) +
  NULL

ggsave_halfheight(
  filename = "plots/hiv-example/p12-only-melding-dists.pdf",
  plot = p3
)
