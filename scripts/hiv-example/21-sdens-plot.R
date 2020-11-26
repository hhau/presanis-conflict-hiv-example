source("scripts/common/plot-settings.R")

library(tibble)
library(magrittr)
library(colorspace)
library(reshape2)
library(wsre)
library(dplyr)

wsre_estimate <- readRDS("rds/hiv-example/big-sub-prior-wsre-est.rds")
phi_grid <- seq(from = 0.01, to = 0.99, length.out = 1000)
n_ests <- length(wsre_estimate$estimates)
all_ests <- wsre_estimate$estimates

res_emp <- lapply(1 : n_ests, function(est_index) {
  tibble(
    phi = phi_grid,
    val = sapply(phi_grid, function(x) sqrt(all_ests[[est_index]]$weighting(x, x))),
    grp = est_index,
    target = "empirical"
  )  
}) %>% bind_rows()

plot_df <- res_emp %>%
  mutate(grp = as.factor(grp))

p1 <- ggplot(plot_df, aes(x = phi, y = val, grp = grp, col = grp)) +
  geom_line() +
  scale_color_manual(
    values = c(qualitative_hcl(20), "#000000")
  ) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

mode_df <- plot_df %>% 
  group_by(grp) %>% 
  transmute(dens_mode = phi[which.max(val)]) %>% 
  distinct() %>% 
  ungroup()

mean_vec <- sapply(all_ests, function(x) x$properties$sample_mean)
lower_vec <- sapply(all_ests, function(x) x$properties$sample_lower_quantile)
upper_vec <- sapply(all_ests, function(x) x$properties$sample_upper_quantile)

mode_df %<>%
  mutate(mean = mean_vec, upper = upper_vec, lower = lower_vec) %>% 
  melt()

vline_df <- mode_df %>% 
  filter(variable %in% c("lower", "upper"))

p2 <- p1 + geom_vline(
  data = vline_df, 
  mapping = aes(xintercept = value, col = grp, lty = variable), 
  alpha = 0.4
) 

# p2
  
ggsave_halfheight(
  filename = "plots/hiv-example/sdens-kde-all.pdf",
  plot = p2
)

## compare to real targets
ref_est <- function(x) dbeta(x, 3.4520971, 0.8341708)

wf_mean <- array(c(wf_mean = seq(from = 0.05, to = 0.8, length.out = 7)))
wf_pars <- list(wf_sd = array(c(wf_sd = 0.0666667 * 1.2)), wf_exponent = 1, target_dimension = 1)

f_gen <- function(a_wf_mu) {
  q_target <- function(x) {
    ref_est(x) * dnorm(x, mean = a_wf_mu, sd = wf_pars$wf_sd)    
  }
  
  nc <- integrate(q_target, 0, 1)
  
  f_target <- function(x) {
    q_target(x) / nc$value
  }
}

res_ana <- lapply(1 : (n_ests - 1), function(est_index) {
  local_wf_mu <- wf_mean[est_index]
  local_target <- f_gen(local_wf_mu)
  res <- tibble(
    phi = phi_grid,
    val = local_target(phi_grid),
    grp = est_index,
    target = "analytic-ish"
  )
}) %>% 
  bind_rows() 

both_df <- bind_rows(res_emp, res_ana) %>% 
  mutate(grp = as.factor(grp), target = as.factor(target)) %>% 
  filter(grp != 8)

p3 <- ggplot(both_df, aes(x = phi, y = val, col = grp, grp = grp, lty = target)) +
  geom_line() + 
  geom_vline(
    data = vline_df %>% filter(grp != 8), 
    mapping = aes(xintercept = value, col = grp),
    lty = 'dashed',
    alpha = 0.4
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave_halfheight(
  filename = "plots/hiv-example/sdens-kde-analytic-compare.pdf",
  plot = p3
)
