source("scripts/common/plot-settings.R")

library(ggplot2)
library(dplyr)

stage_one_samples <- readRDS("rds/hiv-example/stage-one-samples.rds")

trace_data <- mcmc_trace_data(stage_one_samples[, ,c("ww", "p[12]")])

p1 <- ggplot(trace_data %>% filter(parameter == "p[12]"), aes(x = iteration, y = value, col = chain)) +
  geom_line() +
  facet_wrap(vars(parameter))

ggsave_halfheight(
  filename = "plots/hiv-example/stage-one-naive-trace.png",
  plot = p1
)
