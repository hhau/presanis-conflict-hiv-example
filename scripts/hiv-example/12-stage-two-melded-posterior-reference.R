source("scripts/hiv-example/11-stage-two-priors.R")

library(parallel)
library(futile.logger)
library(abind)
library(dplyr)

flog.info("Reading data and stage one samples")

data <- readRDS("rds/hiv-example/hiv-data.rds")
pars_of_interest <- readRDS(
  file = "rds/hiv-example/pars-of-interest.rds"
)

stage_one_samples <- readRDS("rds/hiv-example/stage-one-reference-samples.rds")
stage_one_phi_samples <- as.vector(stage_one_samples[, , "p[12]"])

# mcmc setup
n_mcmc <- 2e4
n_chain <- 5

# stage two log posterior
# stage two prior marginal cancels? whatever I choose for it
# it factors out through the dgp
log_posterior <- function(phi_nu, phi_de) {
  res <-
    dbinom(x = data$y[12], size = data$n[12], prob = phi_nu, log = TRUE) -
    dbinom(x = data$y[12], size = data$n[12], prob = phi_de, log = TRUE)
  return(res)
}

# run through the original samples (non-wsre) and evaluate them under the second
# model.
flog.info("Starting MCMC loop")
log_filename <- paste0("logs/" , Sys.Date(), "-stage-two-run.log") 
flog.appender(appender.file(log_filename), name = "stage-two-logger")

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
  
  for (mcmc_index in 2 : (n_mcmc + 1)) {
    # propose an index,
    index_prop <- sample.int(
      length(stage_one_phi_samples),
      size = 1
    )
    index_curr <- stage_one_indices[mcmc_index - 1, 1, 1]
    # draw phi sample
    phi_prop <- stage_one_phi_samples[index_prop]
    phi_curr <- phi_samples[mcmc_index - 1, 1, 1]

    # calc accept probability,
    log_alpha <- 
      log_posterior(phi_prop, phi_curr) +
      log_logarithmic_pooled_prior_beta_approx(phi_prop, phi_curr)

    # accept/reject
    if (runif(1) < exp(log_alpha)) {
      stage_one_indices[mcmc_index, 1, 1] <- index_prop
      phi_samples[mcmc_index, 1, 1] <- phi_prop
    } else {
      stage_one_indices[mcmc_index, 1, 1] <- index_curr
      phi_samples[mcmc_index, 1, 1] <- phi_curr
    }

    # append to log file.
    if (mcmc_index %% 50 == 0) {
      flog.info(
        sprintf("Chain: %d, Iteration: %d", chain_index, mcmc_index),
        name = "stage-two-logger"
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
  file = "rds/hiv-example/stage-two-reference-samples.rds"
)
