source("scripts/common/plot-settings.R")

library(tibble)
library(dplyr)

wsre_samples <- readRDS(
  file = "rds/hiv-example/stage-two-wsre-samples.rds"
)

no_wsre_samples <- readRDS(
  file = "rds/hiv-example/stage-two-samples.rds"
)

reference_samples <- readRDS(
  file = "rds/hiv-example/stage-two-reference-samples.rds"
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

reference_quantiles <- quantile(
  x = as.vector(reference_samples$phi_samples),
  probs = quantiles_of_interest
)

plot_tbl <- tibble(
  x = rep(reference_quantiles, times = 2),
  y = c(wsre_quantiles, no_wsre_quantiles),
  grp = factor(c(
    rep("B_WSRE", times = length(quantiles_of_interest)),
    rep("A_No-WSRE", times = length(quantiles_of_interest))
  ), ordered = TRUE)
)

p1 <- ggplot(plot_tbl, aes(x = x, y = y, pch = grp)) +
  geom_point(alpha = 1) +
  geom_abline(slope = 1, intercept = 0) +
  labs(
    col = "Approach"
  ) +
  xlab("Reference quantiles") +
  ylab("Melded model quantiles") +
  scale_discrete_manual(
    aesthetics = "pch",
    values = c(
      "B_WSRE" = 1,
      "A_No-WSRE" = 3
    ) ,
    labels = c(
      "B_WSRE" = "WSRE",
      "A_No-WSRE" = "Naive"
    )
  ) +
  scale_x_continuous(
    limits = c(0.15, 0.45)
  ) +
  scale_y_continuous(
    limits = c(0.15, 0.45)
  ) +
  coord_fixed() +
  labs(pch = "Method")

ggsave(
  filename = "plots/hiv-example/posterior-qq-plot.pdf",
  plot = p1,
  width = 12.5,
  height = 10.5,
  units = 'cm'
)

# # ks.testing
# ks.test(
#   x = as.numeric(wsre_samples$phi_samples),
#   y = as.numeric(reference_samples)
# )
# 
# ks.test(
#   x = as.numeric(no_wsre_samples$phi_samples),
#   y = as.numeric(reference_samples)
# )
# 
# ks.test(
#   x = as.numeric(no_wsre_samples$phi_samples),
#   y = as.numeric(wsre_samples$phi_samples)
# )

