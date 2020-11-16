source("scripts/common/plot-settings.R")

plot_df <- readRDS("rds/hiv-example/ratio-estimates-df.rds")
wsre_estimate <- readRDS("rds/hiv-example//big-sub-prior-wsre-est.rds")
naive_prior_samples <- readRDS("rds/hiv-example/prior-samples.rds") 
ref_melded_posterior_samples <- readRDS("rds/hiv-example/stage-two-reference-samples.rds")

p_1 <- ggplot(plot_df, aes(x = x_nu, y = log1p(log1p(value)), col = d_type)) +
  geom_line() +
  facet_wrap(vars(x_de_plot), labeller = label_parsed) +
  labs(col = "Estimate") +
  ggtitle(
    sprintf(
      "N# WSRE: %d, N# Naive: %d", 
      wsre_estimate$properties$n_mcmc_samples, 
      prod(dim(naive_prior_samples)[1 : 2])
    )
  )

ggsave_halfheight(
  filename = "plots/hiv-example/ratio-estimates-full.pdf",
  plot = p_1
)

p_2 <- p_1 + 
  scale_x_continuous(
    limits = c(
      min(ref_melded_posterior_samples$phi_samples),
      max(ref_melded_posterior_samples$phi_samples)
    )
  ) +
  scale_y_sqrt(limits = c(0, 2)) 

ggsave_halfheight(
  filename = "plots/hiv-example/ratio-estimates-zoomed.pdf",
  plot = p_2
)