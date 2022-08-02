source("scripts/common/plot-settings.R")

library(ggplot2)
library(dplyr)
library(colorspace)

stage_one_samples <- readRDS("rds/hiv-example/stage-one-samples.rds")
thin_vec <- round(seq(from = 0, to = dim(stage_one_samples)[1], length.out = 2e3))
thin_join <- tibble(iteration_old = thin_vec, iteration_new = 1 : 2e3)
trace_data <- mcmc_trace_data(stage_one_samples[, ,c("ww", "p[12]")]) %>%
  filter(parameter == "p[12]", iteration %in% thin_vec) %>%
  left_join(thin_join, by = c('iteration' = 'iteration_old')) %>%
  mutate(plot_parameter = "phi", stage = "Stage one") %>%
  select(-iteration) %>%
  rename(iteration = iteration_new)

stage_two_samples <- readRDS("rds/hiv-example/stage-two-samples.rds")
trace_data_two <- mcmc_trace_data(stage_two_samples$phi_samples) %>%
  mutate(plot_parameter = "phi", stage = "Stage two")

combo_df <- bind_rows(trace_data, trace_data_two)

fancy_scientific <- function(l) {
 l <- format(l, scientific = TRUE)
 l <- gsub("1e\\+00","1",l)
 l <- gsub("^(.*)e", "'\\1'e", l)
 l <- gsub("e", "%*%10^", l)
 l <- gsub("\\+", "", l)
 return(parse(text = l))
}

p1 <- ggplot(combo_df, aes(x = iteration, y = value, col = chain)) +
  geom_line() +
  facet_wrap(vars(stage), scales = "fixed", ncol = 1, nrow = 2) +
  scale_color_manual(
    values = c(diverging_hcl(trace_data$n_chains[1]))
  ) +
  labs(col = "Chain") +
  xlab("Iteration") +
  ylab(expression(phi)) +
  # scale_x_log10(
  #   breaks = 10^(0 : 5),
  #   labels = sapply(c("1", sprintf("10^%d", 1 : 5)), function(x) parse(text = x))
  # ) +
  scale_x_continuous(expand = c(0.025, 0.025)) +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    # strip.text = element_text(size = rel(1.5)),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave_base(
  filename = "plots/hiv-example/stage-one-naive-trace.png",
  plot = p1,
  width = 12.5 * 0.72 - 2,
  height = display_settings$half_page_plot_width * 1.5
)
