# checking out the new prior on p_{12}
library(dplyr)

n_sim <- 5e3

tar_beta <- 9

a_samples <- rbeta(n = n_sim, shape1 = 1, shape2 = 2)
b_samples <- rbeta(n = n_sim, shape1 = 1, shape2 = tar_beta)
c_samples <- rbeta(n = n_sim, shape1 = 1, shape2 = tar_beta)
d_samples <- rbeta(n = n_sim, shape1 = 1, shape2 = tar_beta)
e_samples <- rbeta(n = n_sim, shape1 = 1, shape2 = 2)
w_samples <- rbeta(n = n_sim, shape1 = 3, shape2 = 1)

e1ab <- e_samples * (1 - a_samples - b_samples)

p12_samples <- (d_samples * b_samples + w_samples * e1ab) / (d_samples * b_samples + e1ab)
hist(p12_samples[between(p12_samples, 0, 1)], xlim = c(0, 1), breaks = 100, freq = FALSE)

