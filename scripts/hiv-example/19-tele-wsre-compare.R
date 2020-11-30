source("scripts/common/plot-settings.R")

library(bayesplot)
library(dplyr)

base_wsre_samples <- readRDS("rds/hiv-example/stage-one-wsre-samples.rds")
tele_wsre_samples <- readRDS("rds/hiv-example/stage-one-wsre-tele-samples.rds")
ref_samples <- readRDS("rds/hiv-example/stage-one-reference-samples.rds")

base_data <- bayesplot:::prepare_mcmc_array(
  base_wsre_samples,
  character(),
  character(),
  list()
) %>% bayesplot:::melt_mcmc.mcmc_array() %>% 
  mutate(wsre_type = "weighted_avg")

tele_data <- bayesplot:::prepare_mcmc_array(
  tele_wsre_samples,
  character(),
  character(),
  list()
) %>% bayesplot:::melt_mcmc.mcmc_array() %>% 
  mutate(wsre_type = "telescope")

ref_data <- bayesplot:::prepare_mcmc_array(
  ref_samples[, , 1 : 9],
  character(),
  character(),
  list()
) %>% bayesplot:::melt_mcmc.mcmc_array() %>% 
  mutate(wsre_type = "reference")  

plot_data <- bind_rows(base_data, tele_data, ref_data) %>% 
  mutate(
    Chain = as.factor(Chain),
    Parameter = as.factor(Parameter),
    wsre_type = as.factor(wsre_type)
  )

p1 <- ggplot(plot_data) +
  geom_density(aes(x = Value, colour = wsre_type)) + 
  facet_wrap(vars(Parameter), scales = "free")

ggsave_halfheight(
  filename = "plots/hiv-example/tele-vs-avg.pdf",
  plot = p1
)
