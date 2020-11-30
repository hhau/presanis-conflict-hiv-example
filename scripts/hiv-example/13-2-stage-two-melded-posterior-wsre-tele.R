source("scripts/hiv-example/11-stage-two-priors.R")
source("scripts/hiv-example/18-telescoping-tests.R")

library(parallel)
library(futile.logger)
library(abind)
library(dplyr)

flog.info("Reading data and stage one samples")

data <- readRDS("rds/hiv-example/hiv-data.rds")
pars_of_interest <- readRDS(
  file = "rds/hiv-example/pars-of-interest.rds"
)

stage_one_samples <- readRDS("rds/hiv-example/stage-one-wsre-tele-samples.rds")
res <- apply(stage_one_samples, 1:2, function(x) {
  db <- x["dd"] * x["bb"]
  e1ab <- x["ee"] * (1 - x["aa"] - x["bb"])
  p12_prop <- (db + x["ww"] * e1ab) / (db + e1ab)
  return(p12_prop)
})
n_iter <- dim(stage_one_samples)[1]
n_chain <- dim(stage_one_samples)[2]
res <- array(as.vector(res), dim = c(n_iter, n_chain, 1))
old_dimnames <- dimnames(stage_one_samples)
old_dimnames[[3]] <- c(old_dimnames[[3]], "p[12]")
stage_one_samples <- abind(stage_one_samples, res, along = 3)
dimnames(stage_one_samples) <- old_dimnames

stage_one_phi_samples <- as.vector(stage_one_samples[, , "p[12]"])

# mcmc setup
n_mcmc <- 1e3
n_chain <- 5

# stage two log posterior
# stage two prior marginal cancels? whatever I choose for it
# it factors out through the dgp
# doesn't really matter, not of wsre interest
log_posterior <- function(phi_nu, phi_de) {
  res <-
    dbinom(x = data$y[12], size = data$n[12], prob = phi_nu, log = TRUE) -
    dbinom(x = data$y[12], size = data$n[12], prob = phi_de, log = TRUE)
  return(res)
}

# run through the original samples (non-wsre) and evaluate them under the second
# model.
flog.info("Starting MCMC loop")
log_filename <- paste0("logs/" , Sys.Date(), "-stage-two-wsre-tele-run.log") 
flog.appender(appender.file(log_filename), name = "stage-two-tele-logger")

mcmc_res <- mclapply(1 : n_chain, mc.cores = n_chain, function(chain_index) {
  # setup containers and initial values
  stage_one_indices <- array(
    data = NA,
    dim = c(n_mcmc + 1, 1, 1),
    dimnames = list(
      sprintf("iteration_%d", 1 : (n_mcmc + 1)),
      sprintf("Chain_%d", chain_index),
      c("phi_1")
    )
  )
  
  phi_samples <- array(
    NA,
    dim = c(n_mcmc + 1, 1, 1),
    dimnames = list(
      sprintf("iteration_%d", 1 : (n_mcmc + 1)),
      sprintf("Chain_%d", chain_index),
      c("phi_1")
    )
  )
  
  stage_one_indices[1, 1, 1] <- sample.int(
    length(stage_one_phi_samples),
    size = 1
  )
  
  phi_samples[1, 1, 1] <- stage_one_phi_samples[
    stage_one_indices[1, 1, 1]
  ]
  
  # n_tele_terms_vec <- array(NA, dim = c(n_mcmc))

  for (ii in 2 : (n_mcmc + 1)) {
    # propose an index,
    index_prop <- sample.int(
      length(stage_one_phi_samples),
      size = 1
    )
    index_curr <- stage_one_indices[ii - 1, 1, 1]
    # draw phi sample
    phi_prop <- stage_one_phi_samples[index_prop]
    phi_curr <- phi_samples[ii - 1, 1, 1]

    # calc accept probability,
    log_alpha <- 
      log_posterior(phi_prop, phi_curr) +
      log_logarithmic_pooled_prior_with_wsre_tele(phi_curr, phi_prop)

    # accept/reject
    if (runif(1) < exp(log_alpha)) {
      stage_one_indices[ii, 1, 1] <- index_prop
      phi_samples[ii, 1, 1] <- phi_prop
    } else {
      stage_one_indices[ii, 1, 1] <- index_curr
      phi_samples[ii, 1, 1] <- phi_curr
    }

    # append to log file.
    if (ii %% 50 == 0) {
      flog.info(
        sprintf("Chain: %d, Iteration: %d", chain_index, ii),
        name = "stage-two-tele-logger"
      )
    }
  }
  return(list(
    phi_samples = phi_samples,
    stage_one_indices = stage_one_indices
  ))
})

flog.info("Finished MCMC loop")

list_subnames <- c("phi_samples", "stage_one_indices")
mcmc_final_res <- lapply(list_subnames, function(a_name) {
  lapply(mcmc_res, function(a_chain_sublist) {
    a_chain_sublist[[a_name]]
  }) %>% 
    abind(along = 2)
})
names(mcmc_final_res) <- list_subnames

flog.info("Writing samples to disk")

saveRDS(
  object = mcmc_final_res,
  file = "rds/hiv-example/stage-two-wsre-tele-samples.rds"
)

# plot_df <- table(n_tele_terms_vec) %>% as.data.frame()
# p1 <- ggplot(plot_df, aes(x = n_tele_terms_vec, y = Freq)) +
#   geom_col()
# 
# ggsave_halfheight(
#   filename = "plots/hiv-example/n-tele-terms-stage-two.pdf",
#   plot = p1
# )
