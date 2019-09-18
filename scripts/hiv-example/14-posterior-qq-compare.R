source("scripts/common/plot-settings.R")

library(tibble)
library(dplyr)

wsre_samples <- readRDS(
  file = "rds/hiv-example/stage-two-wsre-samples.rds"
)

no_wsre_samples <- readRDS(
  file = "rds/hiv-example/stage-two-samples.rds"
)

joint_samples <- readRDS(
  file = "rds/hiv-example/full-model-fit.rds"
)

quantiles_of_interest <- seq(from = 0.005, to = 0.995, by = 0.005)

wsre_quantiles <- quantile(
  x = as.vector(wsre_samples$phi_samples),
  probs = quantiles_of_interest
)

no_wsre_quantiles <- quantile(
  x = as.vector(no_wsre_samples$phi_samples),
  probs = quantiles_of_interest
)

joint_quantiles <- quantile(
  x = as.vector(joint_samples[ , , "p[12]"]),
  probs = quantiles_of_interest
)

plot_tbl <- tibble(
  x = rep(joint_quantiles, times = 2),
  y = c(wsre_quantiles, no_wsre_quantiles),
  grp = factor(c(
    rep("a_WSRE", times = length(quantiles_of_interest)),
    rep("b_No-WSRE", times = length(quantiles_of_interest))
  ), ordered = TRUE)
)

p1 <- ggplot(plot_tbl, aes(x = x, y = y, col = grp)) +
  geom_point(alpha = 1) +
  geom_abline(slope = 1, intercept = 0) +
  labs(
    col = "Method"
  ) +
  xlab("Joint model quantiles") +
  ylab("Melded model quantiles") +
  scale_discrete_manual(
    aesthetics = "col",
    values = c(
      "a_WSRE" = as.character(blues[2]),
      "b_No-WSRE" = highlight_col
    ) ,
    labels = c(
      "a_WSRE" = "WSRE",
      "b_No-WSRE" = "No WSRE"
    )
  )

ggsave_halfheight(
  filename = "plots/hiv-example/posterior-qq-plot.pdf",
  plot = p1
)

