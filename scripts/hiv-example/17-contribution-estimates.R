source("scripts/common/plot-settings.R")

library(wsre)
library(dplyr)
library(colorspace)
library(ggplot2)
library(wsre)
library(tibble)

wsre_estimate <- readRDS("rds/hiv-example/big-sub-prior-wsre-est.rds")

# code to make proportion plots at same x_de values, for a small grid of 
# x_nu values
x_de_values <- seq(from = 0.15, to = 0.75, length.out = 6)
x_nu_grid_values <- seq(from = 0.05, to = 0.95, length.out = 10)

x_input <- expand.grid(x_nu = x_nu_grid_values, x_de = x_de_values)

temp_func <- function(a_nu_value, a_de_value) {
  n_estimates <- length(wsre_estimate$estimates)
  res <- list()
  for (ii in 1 : n_estimates) {
    an_est <- wsre_estimate$estimates[[ii]]
    res[[ii]] <- tibble(
      x_nu = a_nu_value,
      x_de = a_de_value,
      wf_index = ii,
      wf_mean = an_est$wf_pars$wf_mean,
      wf_value = an_est$weighting(a_nu_value, a_de_value),
      sample_mean = an_est$properties$sample_mean,
      wf_bw = an_est$properties$bandwidth
    )
  }
  r_res <- bind_rows(res)
  r_res$normalised_proportion <- r_res$wf_value / sum(r_res$wf_value)
  return(r_res)
}

plot_df <- lapply(1 : nrow(x_input), function(ii) {
  temp_func(
    x_input[ii, ]$x_nu,
    x_input[ii, ]$x_de
  )
}) %>% 
  bind_rows()

estimate_label_makes <- function(an_est) {
  lab_str <- sprintf(
    "list(mu['w']: ~ %s, ~ bar(phi): ~ %s, ~ hat(sigma)[phi]: ~ %s, ~ italic('h'): ~ %s)",
    an_est$wf_pars$wf_mean %>% format(digits = 3, nsmall = 3),
    an_est$properties$sample_mean %>% format(digits = 3, nsmall = 3),
    an_est$properties$sample_sd %>% format(digits = 3, nsmall = 3),
    an_est$properties$bandwidth %>% format(digits = 3, nsmall = 4) 
  ) %>% as.expression()
}

recode_vec <- sapply(wsre_estimate$estimates, estimate_label_makes)

plot_df <- plot_df %>% 
  mutate(
    wf_index = factor(
      wf_index,
      levels = unique(wf_index),
      labels = recode_vec
    ), 
    x_nu = as.factor(x_nu)
  )

plot_df$x_de_plot <- sprintf("'x'['de'] ~ '=' ~ %.2f", plot_df$x_de)

p_1 <- ggplot(plot_df, aes(x = x_nu, fill = wf_index)) +
  facet_wrap(vars(x_de_plot), labeller = label_parsed) +
  geom_col(aes(y = normalised_proportion), colour = "white") +
  theme(
    axis.text.x = element_text(angle = 90)
  ) + 
    scale_fill_manual(
      labels = function(x) parse(text = x),
      values = c(qualitative_hcl(12), "#000000")
    ) +
  labs(fill = expression(mu['w']))

p_1

ggsave_base(
  filename = "plots/hiv-example/contributing-wf.pdf",
  plot = p_1, 
  width = 21,
  height = 15
)
