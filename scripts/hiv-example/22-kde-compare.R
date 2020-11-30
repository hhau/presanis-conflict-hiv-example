library(wsre)
library(dplyr)
library(parallel)
library(pbapply)
library(tibble)
library(patchwork)

ref_samples <- readRDS("rds/hiv-example/reference-prior-samples.rds")
ref_phi_samples <- ref_samples[, , "p[12]"] %>% as.numeric()
sample_sizes <- seq(from = 500, to = 5000, by = 500)
x_vals <- seq(from = 0, to = 1, length.out = 250)
quantiles_of_interest <- c(0.05, 0.1, 0.5, 0.9, 0.95)

n_mc <- 100
input_grid <- expand.grid(mc_index = 1 : n_mc, sample_size = sample_sizes)

res <- pblapply(X = 1 : nrow(input_grid), cl = 6, FUN = function(row_index) {
  local_sample_size <- input_grid[row_index, "sample_size"]
  local_mc_index <- input_grid[row_index, "mc_index"]
  
  phi_prior_samples <- sample(
    x = ref_phi_samples,
    size = local_sample_size,
    replace = TRUE
  ) %>% as.matrix()
  
  kde_bw <- bw.SJ(phi_prior_samples)
  curve_f <- function(x) {
    log(wsre:::kde_func_nd(x, phi_prior_samples, kde_bw))
  }
  
  plot_df <- data.frame(
    x = x_vals,
    y = sapply(x_vals, curve_f),
    sample_size = local_sample_size,
    mc_index = local_mc_index
  )  
})

plot_df <- bind_rows(res) %>% 
  as_tibble()

p1 <- ggplot(plot_df, aes(x = x, y = y, group = mc_index)) +
  geom_line(alpha = 0.1) +
  facet_wrap(vars(sample_size), scales = "free_x") +
  scale_x_continuous(limits = c(0, 0.15)) +
  theme(axis.title.y = element_blank())

uncert_df <- plot_df %>% 
  group_by(x, sample_size) %>% 
  mutate(
    median = quantile(y, 0.5),
    q10 = quantile(y, 0.1),
    q90 = quantile(y, 0.9),
    sample_size = as.factor(sample_size)
  )

p2 <- ggplot(uncert_df, aes(x = x, group = sample_size)) +
  geom_line(aes(y = median)) +
  geom_ribbon(aes(ymin = q10, ymax = q90), alpha = 0.2) +
  facet_wrap(vars(sample_size), scales = "free_x") +
  scale_x_continuous(limits = c(0, 0.15)) +
  theme(axis.title.y = element_blank())

p1 + p2
