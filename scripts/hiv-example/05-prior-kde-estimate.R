library(dplyr)
library(wsre)

prior_samples <- readRDS(file = "rds/hiv-example/prior-samples.rds")
index_of_interest <- grep("12", dimnames(prior_samples)[["parameters"]])
p12_samples <- prior_samples[, , index_of_interest] %>% as.vector()
p12_bw <- bw.SJ(p12_samples)

log_prior_density_est <- function(x) {
  log(wsre:::gauss_kde(x, p12_samples, p12_bw))
}

curve(log_prior_density_est)
