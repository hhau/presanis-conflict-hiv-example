source("scripts/common/plot-settings.R")

library(dplyr)

plot_df <- readRDS("rds/hiv-example/ratio-estimates-df.rds")
wsre_estimate <- readRDS("rds/hiv-example/big-sub-prior-wsre-est.rds")
naive_prior_samples <- readRDS("rds/hiv-example/prior-samples.rds") 
ref_melded_posterior_samples <- readRDS("rds/hiv-example/stage-two-reference-samples.rds")

n_ests <- length(wsre_estimate$estimates)
all_ests <- wsre_estimate$estimates
mean_vec <- sapply(all_ests, function(x) x$properties$sample_mean) %>% data.frame(xintercept = .)

p_1 <- ggplot(plot_df, aes(x = x_nu, y = log1p(log1p(value)), col = d_type)) +
  geom_line() +
  facet_wrap(vars(x_de_plot), labeller = label_parsed) +
  labs(col = "Estimate") +
  ggtitle(
    sprintf(
      "N# WSRE: %0.2f, N# Naive: %d", 
      wsre_estimate$properties$n_mcmc_samples, 
      prod(dim(naive_prior_samples)[1 : 2])
    )
  ) +
  geom_vline(data = mean_vec, aes(xintercept = xintercept), alpha = 0.2, lty = "dashed") +
  xlab(expression(x['nu']))

ggsave_base(
  filename = "plots/hiv-example/ratio-estimates-full.pdf",
  plot = p_1,
  width = 21,
  height = 15
)

p_2 <- p_1 + 
  scale_x_continuous(
    limits = c(
      min(ref_melded_posterior_samples$phi_samples),
      max(ref_melded_posterior_samples$phi_samples)
    )
  ) +
  scale_y_sqrt(limits = c(0, 2)) 

ggsave_base(
  filename = "plots/hiv-example/ratio-estimates-zoomed.pdf",
  plot = p_2,
  width = 21,
  height = 15
)

inverse_plot_df <- readRDS("rds/hiv-example/ratio-estimates-df-inverse.rds")

p_3 <- ggplot(inverse_plot_df, aes(x = x_de, y = log1p(log1p(value)), col = d_type)) +
  geom_line() +
  facet_wrap(vars(x_nu_plot), labeller = label_parsed) +
  labs(col = "Estimate") +
  ggtitle(
    sprintf(
      "N# WSRE: %0.2f, N# Naive: %d", 
      wsre_estimate$properties$n_mcmc_samples, 
      prod(dim(naive_prior_samples)[1 : 2])
    )
  ) +
  geom_vline(data = mean_vec, aes(xintercept = xintercept), alpha = 0.2, lty = "dashed") +
  xlab(expression(x['de']))

# inverse_plot_df %>% 
#   filter(d_type %in% c("ref", "naive", "wsre_tele_fixed_N", "wsre_tele1d_no_naive")) %>% 
# ggplot(aes(x = x_de, y = log1p(log1p(value)), col = d_type)) +
#   geom_line() +
#   facet_wrap(vars(x_nu_plot), labeller = label_parsed) +
#   labs(col = "Estimate") +
#   ggtitle(
#     sprintf(
#       "N# WSRE: %0.2f, N# Naive: %d", 
#       wsre_estimate$properties$n_mcmc_samples, 
#       prod(dim(naive_prior_samples)[1 : 2])
#     )
#   ) +
#   geom_vline(data = mean_vec, aes(xintercept = xintercept), alpha = 0.2, lty = "dashed") +
#   xlab(expression(x['de']))

ggsave_base(
  filename = "plots/hiv-example/ratio-estimates-full-inverse.pdf",
  plot = p_3,
  width = 21,
  height = 15
)
