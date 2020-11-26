library(tibble)
library(dplyr)
library(wsre)
library(pbapply)
library(magrittr)

source("scripts/hiv-example/18-telescoping-tests.R")
source("scripts/common/plot-settings.R")

ref_stage_one_samples <- readRDS(
  "rds/hiv-example/stage-one-reference-samples.rds"
)

ref_melded_posterior_samples <- readRDS(
  "rds/hiv-example/stage-two-reference-samples.rds"
)

wsre_estimate <- readRDS(
  "rds/hiv-example/big-sub-prior-wsre-est.rds"
)

naive_prior_samples <- readRDS(
  file = "rds/hiv-example/prior-samples.rds"
)

naive_kde_bw <- bw.SJ(naive_prior_samples[, , "p[12]"] %>% as.numeric())
naive_kde <- function(x) {
  wsre:::kde_func_nd(
    x, 
    naive_prior_samples[, , "p[12]"] %>% 
      as.numeric() %>% 
      as.matrix()
    , naive_kde_bw
  )
}

naive_ratio <- function(x_nu, x_de) {
  naive_kde(x_nu) / naive_kde(x_de)
}

ref_ratio <- function(x_nu, x_de) {
  dbeta(x_nu, 3.4520971, 0.8341708) / dbeta(x_de, 3.4520971, 0.8341708) 
}

# under uniform values for x_nu and x_de
n_points_to_test <- 5000

all_x_nu_unif <- runif(n_points_to_test)
all_x_de_unif <- runif(n_points_to_test)

naive_errors_unif <- array(NA, dim = n_points_to_test)
wavg_errors_unif <- array(NA, dim = n_points_to_test)
tele_errors_unif <- array(NA, dim = n_points_to_test)
tele_dumb_errors_unif <- array(NA, dim = n_points_to_test)

for (ii in 1 : n_points_to_test) {
  ref_val <- ref_ratio(all_x_nu_unif[ii], all_x_de_unif[ii])
  naive_val <- naive_ratio(all_x_nu_unif[ii], all_x_de_unif[ii])
  wavg_val <- evaluate(wsre_estimate, all_x_nu_unif[ii], all_x_de_unif[ii])
  tele_val <- evaluate_telescope_1d(wsre_estimate, all_x_nu_unif[ii], all_x_de_unif[ii])
  tele_dumb_val <- evaluate_telescope_1d_dumb(wsre_estimate, all_x_nu_unif[ii], all_x_de_unif[ii])
  
  naive_errors_unif[ii] <- (naive_val - ref_val)^2
  wavg_errors_unif[ii] <- (wavg_val - ref_val)^2
  tele_errors_unif[ii] <- (tele_val - ref_val)^2
  tele_dumb_errors_unif[ii] <- t(tele_dumb_val - ref_val)^2
}

plot_df <- tibble(
  x = c(naive_errors_unif, wavg_errors_unif, tele_errors_unif, tele_dumb_errors_unif),
  est_type = rep(c("naive", "wavg", "tele", "tele_dumb"), each = n_points_to_test),
  point_dist = rep("unif", n_points_to_test * 4)
)

# stage 1 values
all_x_nu_stage_1 <- ref_stage_one_samples[, , "p[12]"] %>% 
  as.numeric() %>% 
  sample(size = n_points_to_test, replace = TRUE)

all_x_de_stage_1 <- ref_stage_one_samples[, , "p[12]"] %>% 
  as.numeric() %>% 
  sample(size = n_points_to_test, replace = TRUE)

naive_errors_stage_1 <- array(NA, dim = n_points_to_test)
wavg_errors_stage_1 <- array(NA, dim = n_points_to_test)
tele_errors_stage_1 <- array(NA, dim = n_points_to_test)
tele_dumb_errors_stage_1 <- array(NA, dim = n_points_to_test)

for (ii in 1 : n_points_to_test) {
  ref_val <- ref_ratio(all_x_nu_stage_1[ii], all_x_de_stage_1[ii])
  naive_val <- naive_ratio(all_x_nu_stage_1[ii], all_x_de_stage_1[ii])
  wavg_val <- evaluate(wsre_estimate, all_x_nu_stage_1[ii], all_x_de_stage_1[ii])
  tele_val <- evaluate_telescope_1d(
    wsre_estimate,
    all_x_nu_stage_1[ii], 
    all_x_de_stage_1[ii]
  )
  tele_dumb_val <- evaluate_telescope_1d_dumb(
    wsre_estimate,
    all_x_nu_stage_1[ii], 
    all_x_de_stage_1[ii]
  )
  
  naive_errors_stage_1[ii] <- (naive_val - ref_val)^2
  wavg_errors_stage_1[ii] <- (wavg_val - ref_val)^2
  tele_errors_stage_1[ii] <- (tele_val - ref_val)^2
  tele_dumb_errors_stage_1[ii] <- (tele_dumb_val - ref_val)^2
}

stage_1_df <- tibble(
  x = c(naive_errors_stage_1, wavg_errors_stage_1, tele_errors_stage_1, tele_dumb_errors_stage_1),
  est_type = rep(c("naive", "wavg", "tele", "tele_dumb"), each = n_points_to_test),
  point_dist = rep("stage_1", n_points_to_test * 4)
)

# stage 2 values
all_x_nu_stage_2 <- ref_melded_posterior_samples %>%
  extract2("phi_samples") %>% 
  as.numeric() %>% 
  sample(size = n_points_to_test, replace = TRUE)

all_x_de_stage_2 <- ref_melded_posterior_samples %>%
  extract2("phi_samples") %>% 
  as.numeric() %>% 
  sample(size = n_points_to_test, replace = TRUE)

naive_errors_stage_2 <- array(NA, dim = n_points_to_test)
wavg_errors_stage_2 <- array(NA, dim = n_points_to_test)
tele_errors_stage_2 <- array(NA, dim = n_points_to_test)
tele_dumb_errors_stage_2 <- array(NA, dim = n_points_to_test)

for (ii in 1 : n_points_to_test) {
  ref_val <- ref_ratio(all_x_nu_stage_2[ii], all_x_de_stage_2[ii])
  naive_val <- naive_ratio(all_x_nu_stage_2[ii], all_x_de_stage_2[ii])
  wavg_val <- evaluate(wsre_estimate, all_x_nu_stage_2[ii], all_x_de_stage_2[ii])
  tele_val <- evaluate_telescope_1d(
    wsre_estimate,
    all_x_nu_stage_2[ii], 
    all_x_de_stage_2[ii]
  )
  
  tele_dumb_val <- evaluate_telescope_1d_dumb(
    wsre_estimate,
    all_x_nu_stage_2[ii], 
    all_x_de_stage_2[ii]
  )
  naive_errors_stage_2[ii] <- (naive_val - ref_val)^2
  wavg_errors_stage_2[ii] <- (wavg_val - ref_val)^2
  tele_errors_stage_2[ii] <- (tele_val - ref_val)^2
  tele_dumb_errors_stage_2[ii] <- (tele_dumb_val - ref_val)^2
}

stage_2_df <- tibble(
  x = c(naive_errors_stage_2, wavg_errors_stage_2, tele_errors_stage_2, tele_dumb_errors_stage_2),
  est_type = rep(c("naive", "wavg", "tele", "tele_dumb"), each = n_points_to_test),
  point_dist = rep("stage_2", n_points_to_test * 4)
)

final_df <- bind_rows(plot_df, stage_1_df, stage_2_df) %>% 
  filter(x > 1e-15) # some tiny things that are not important

median_df <- final_df %>% 
  group_by(est_type, point_dist) %>% 
  transmute(median = median(x)) %>% 
  distinct()

all_plot <- ggplot(final_df, aes(x = x)) +
  geom_density() +
  geom_vline(
    data = median_df, 
    mapping = aes(xintercept = median), 
    col = "red",
    alpha = 0.5
  ) +
  geom_rug(alpha = 0.2) +
  facet_grid(rows = vars(est_type), cols = vars(point_dist), scales = "free_x") +
  scale_x_log10() +
  xlab("Squared error")

ggsave_fullpage(
  filename = "plots/hiv-example/sq-error-distribution.pdf",
  plot = all_plot
)

sub_df <- final_df %>% 
  filter(point_dist == "stage_2", est_type != "naive") 

sub_vline_df <- median_df %>% 
  filter(point_dist == "stage_2", est_type != "naive")

small_plot <- ggplot(sub_df, aes(x = x, col = est_type)) +
  geom_density() +
  geom_vline(
    data = sub_vline_df, 
    mapping = aes(xintercept = median, col = est_type), 
    alpha = 0.5
  ) +
  scale_x_log10() +
  geom_rug(alpha = 0.2) +
  ggtitle("Distribution of squared errors, evaluated at Stage 2 points.") +
  xlab("Squared error") 
  
ggsave_halfheight(
  filename = "plots/hiv-example/stage-two-wsre-sq-error-distributions.pdf",
  plot = small_plot
)
