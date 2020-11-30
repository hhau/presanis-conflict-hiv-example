source("scripts/common/plot-settings.R")

library(patchwork)
library(tibble)
library(dplyr)
library(colorspace)

stage_one_samples <- readRDS("rds/hiv-example/stage-one-samples.rds")
wsre_samples <- readRDS("rds/hiv-example/stage-two-wsre-samples.rds")
wsre_tele_samples <- readRDS("rds/hiv-example/stage-two-wsre-tele-samples.rds")
no_wsre_samples <- readRDS("rds/hiv-example/stage-two-samples.rds")
reference_samples <- readRDS("rds/hiv-example/stage-two-reference-samples.rds")

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
  x = rep(reference_quantiles, times = 2),
  y = c(wsre_quantiles, no_wsre_quantiles),
  grp = factor(c(
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
  scale_colour_manual(
    name = "Method",
    values = c(
      "B_WSRE" = blues[2],
      "A_No-WSRE" = highlight_col
    ) ,
    labels = c(
      "B_WSRE" = "WSRE",
      "A_No-WSRE" = "Naive"
    )
  ) +
  scale_x_continuous(
    limits = c(0.0, 0.4)
  ) +
  scale_y_continuous(
    limits = c(0.0, 0.4)
  ) +
  coord_fixed() +
  guides(
    pch = guide_legend(reverse = TRUE), 
    col = guide_legend(reverse = TRUE)
  ) 

ggsave(
  filename = "plots/hiv-example/posterior-qq-plot.pdf",
  plot = p1,
  width = 12.5 * 0.72,
  height = 10.5 * 0.72,
  units = 'cm'
)

fancy_scientific <- function(l) {
 l <- format(l, scientific = TRUE)
 l <- gsub("1e\\+00","1",l)
 l <- gsub("^(.*)e", "'\\1'e", l)
 l <- gsub("e", "%*%10^", l)
 l <- gsub("+", "", l)
 return(parse(text = l))
}

trace_data <- mcmc_trace_data(stage_one_samples[, ,c("ww", "p[12]")]) %>% 
  filter(parameter == "p[12]") %>% 
  mutate(plot_parameter = "phi")

ptrace <- ggplot(trace_data, aes(x = iteration, y = value, col = chain)) +
  geom_line() +
  facet_wrap(vars(plot_parameter), labeller = label_parsed) +
  scale_color_manual(
    values = c(diverging_hcl(trace_data$n_chains[1]))
  ) +
  labs(col = "Chain") +
  xlab(expression(log[10]("Iteration"))) + 
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, margin = margin(10, 10, 0, 0)),
    strip.text = element_text(size = rel(1.25)),
    legend.position = "none"
  ) + 
  scale_x_log10(labels = fancy_scientific)

pfinal <- p1 + ptrace
pvert_final <- ptrace / p1

ggsave_halfheight(
  filename = "plots/hiv-example/stage-one-trace-and-posterior-qq-plot.png",
  plot = pfinal
)

ggsave_base(
  filename = "plots/hiv-example/stage-one-trace-and-posterior-qq-plot-vertical.png",
  plot = pvert_final,
  height = display_settings$full_page_plot_width,
  width = display_settings$half_page_plot_width
)
