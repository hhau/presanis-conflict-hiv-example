source("scripts/common/plot-settings.R")

library(ggplot2)
library(dplyr)
library(colorspace)

stage_one_samples <- readRDS("rds/hiv-example/stage-one-samples.rds")
trace_data <- mcmc_trace_data(stage_one_samples[, ,c("ww", "p[12]")]) %>% 
  filter(parameter == "p[12]") %>% 
  mutate(plot_parameter = "phi")

fancy_scientific <- function(l) {
 l <- format(l, scientific = TRUE)
 l <- gsub("1e\\+00","1",l)
 l <- gsub("^(.*)e", "'\\1'e", l)
 l <- gsub("e", "%*%10^", l)
 l <- gsub("+", "", l)
 return(parse(text = l))
}

p1 <- ggplot(trace_data, aes(x = iteration, y = value, col = chain)) +
  geom_line() +
  facet_wrap(vars(plot_parameter), labeller = label_parsed) +
  scale_color_manual(
    values = c(diverging_hcl(trace_data$n_chains[1]))
  ) +
  labs(col = "Chain") +
  xlab("Iteration") + 
  scale_x_log10(labels = fancy_scientific) +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45, margin = margin(10, 10, 0, 0)),
    strip.text = element_text(size = rel(1.5)),
    legend.position = "none"
  ) 

ggsave_base(
  filename = "plots/hiv-example/stage-one-naive-trace.png",
  plot = p1,
  width = display_settings$half_page_plot_width + 2,
  height = display_settings$half_page_plot_width
)
