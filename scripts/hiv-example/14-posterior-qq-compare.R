source("scripts/common/plot-settings.R")

library(tibble)
library(dplyr)

wsre_samples <- readRDS(
  file = "rds/hiv-example/stage-two-wsre-samples.rds"
)

wsre_tele_samples <- readRDS(
  file = "rds/hiv-example/stage-two-wsre-tele-samples.rds"
)

no_wsre_samples <- readRDS(
  file = "rds/hiv-example/stage-two-samples.rds"
)

reference_samples <- readRDS(
  file = "rds/hiv-example/stage-two-reference-samples.rds"
)

quantiles_of_interest <- seq(from = 0.01, to = 0.99, by = 0.01)

wsre_quantiles <- quantile(
  x = as.vector(wsre_samples$phi_samples),
  probs = quantiles_of_interest
)

wsre_tele_quantiles <- quantile(
  x = as.vector(wsre_tele_samples$phi_samples),
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
  x = rep(reference_quantiles, times = 3),
  y = c(wsre_tele_quantiles, wsre_quantiles, no_wsre_quantiles),
  grp = factor(c(
    rep("C_WSRE_TELE", times = length(quantiles_of_interest)),
    rep("B_WSRE", times = length(quantiles_of_interest)),
    rep("A_No-WSRE", times = length(quantiles_of_interest))
  ), ordered = TRUE)
)

# line version
p1 <- ggplot(plot_tbl, aes(x = x, y = y, col = grp)) +
  geom_line(alpha = 0.75, size = rel(0.5)) +
  geom_abline(slope = 1, intercept = 0) +
  labs(
    col = "Approach"
  ) +
  xlab("Reference quantiles") +
  ylab("Melded model quantiles") +
  # scale_shape_manual(
  #   name = "Method",
  #   values = c(
  #     "C_WSRE_TELE" = 15,
  #     "B_WSRE" = 19,
  #     "A_No-WSRE" = 4
  #   ) ,
  #   labels = c(
  #     "C_WSRE_TELE" = "WSRE (Telescoping)",
  #     "B_WSRE" = "WSRE",
  #     "A_No-WSRE" = "Naive"
  #   )
  # ) +
  scale_colour_manual(
    name = "Method",
    values = c(
      "C_WSRE_TELE" = blues[1],
      "B_WSRE" = blues[2],
      "A_No-WSRE" = highlight_col
    ) ,
    labels = c(
      "C_WSRE_TELE" = "WSRE (Telescoping)",
      "B_WSRE" = "WSRE",
      "A_No-WSRE" = "Naive"
    )
  ) +
  # scale_x_continuous(
  #   limits = c(0.15, 0.4)
  # ) +
  # scale_y_continuous(
  #   limits = c(0.15, 0.4)
  # ) +
  coord_fixed() +
  guides(pch = guide_legend(reverse = TRUE), col = guide_legend(reverse = TRUE))

# point version
# p1 <- ggplot(plot_tbl, aes(x = x, y = y, shape = grp, col = grp)) +
#   geom_point(alpha = 0.85, size = rel(1.25)) +
#   geom_abline(slope = 1, intercept = 0) +
#   labs(
#     col = "Approach"
#   ) +
#   xlab("Reference quantiles") +
#   ylab("Melded model quantiles") +
#   scale_shape_manual(
#     name = "Method",
#     values = c(
#       "C_WSRE_TELE" = 15,
#       "B_WSRE" = 19,
#       "A_No-WSRE" = 4
#     ) ,
#     labels = c(
#       "C_WSRE_TELE" = "WSRE (Telescoping)",
#       "B_WSRE" = "WSRE",
#       "A_No-WSRE" = "Naive"
#     )
#   ) +
#   scale_colour_manual(
#     name = "Method",
#     values = c(
#       "C_WSRE_TELE" = blues[1],
#       "B_WSRE" = blues[2],
#       "A_No-WSRE" = highlight_col
#     ) ,
#     labels = c(
#       "C_WSRE_TELE" = "WSRE (Telescoping)",
#       "B_WSRE" = "WSRE",
#       "A_No-WSRE" = "Naive"
#     )
#   ) +
#   scale_x_continuous(
#     limits = c(0.15, 0.4)
#   ) +
#   scale_y_continuous(
#     limits = c(0.15, 0.4)
#   ) +
#   coord_fixed() +
#   guides(pch = guide_legend(reverse = TRUE), col = guide_legend(reverse = TRUE))

ggsave(
  filename = "plots/hiv-example/posterior-qq-plot.pdf",
  plot = p1,
  width = 12.5,
  height = 10.5,
  units = 'cm'
)
