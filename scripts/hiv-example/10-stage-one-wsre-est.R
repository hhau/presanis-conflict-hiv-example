source("scripts/common/plot-settings.R")

library(wsre)
library(rstan)
library(parallel)
library(futile.logger) # substantial overhead for writing to log file
library(logitnorm)
library(abind)

# load data
data <- readRDS("rds/hiv-example/hiv-data.rds")
pars_of_interest <- readRDS(
  file = "rds/hiv-example/pars-of-interest.rds"
)

# Precomplie stan model / load wsre est
prefit <- stan_model("scripts/stan-files/hiv-ev-sythn-big-submodel.stan")
stan_data <- list(
  n_studies = nrow(data),
  y_obs = data$y,
  n_obs = data$n
)

wsre_est <- readRDS(
  file = "rds/hiv-example/big-sub-prior-wsre-est.rds"
)

# need to do this to generate a stanfit object 
# - which is the one with the log prob method
model_fit <- suppressWarnings(sampling(
  prefit,
  data = stan_data,
  chains = 1,
  iter = 1,
  refresh = 0
))

par_names <- c("aa", "bb", "cc", "dd", "ee", "ff", "gg", "hh", "ww")
n_par <- length(par_names)

# setup stage one containers
n_iter <- (50000 + 2000)
n_chain <- 6
thin_vec <- round(seq(from = 2000, to = n_iter, length.out = 2000))

# tuning parameters
prop_sigma <- c(
  aa = 0.05,
  bb = 0.15,
  cc = 0.05,
  dd = 0.1,
  ee = 0.1,
  ff = 0.25,
  gg = 0.25,
  hh = 0.25,
  ww = 0.125
)

# (try mclapply? I know this occasionally has issues stan):
res_list <- mclapply(1 : n_chain, mc.cores = n_chain, function(chain_id) {
  # set up containers
  # use bayesplot format [iteration, chain, parameter]
  # keep a defunct index to abind over later?
  res_array <- array(NA, dim = c(n_iter, 1, n_par))
  dimnames(res_array) <- list(
    sprintf("iteration-%d", 1 : n_iter),
    sprintf("chain-%d", chain_id),
    par_names
  )
  
  log_file_name <- sprintf("logs/mcmc-stage-one-wsre-log.log")
  flog.appender(appender.file(log_file_name), "wsre-logger")

  # initialise
  res_array[1, 1, ] <- c(
    aa = rlogitnorm(n = 1, mu = logit(0.1), sigma = 0.05),
    bb = rlogitnorm(n = 1, mu = logit(0.012), sigma = 0.05),
    cc = rlogitnorm(n = 1, mu = logit(0.016), sigma = 0.05),
    dd = rlogitnorm(n = 1, mu = logit(0.025), sigma = 0.05),
    ee = rlogitnorm(n = 1, mu = logit(0.00045), sigma = 0.05),
    ff = rlogitnorm(n = 1, mu = logit(0.3), sigma = 0.125),
    gg = rlogitnorm(n = 1, mu = logit(0.5), sigma = 0.25),
    hh = rlogitnorm(n = 1, mu = logit(0.5), sigma = 0.25),
    ww = rlogitnorm(n = 1, mu = logit(0.125), sigma = 0.125)
  )

  # loop
  for (ii in 2 : (n_iter)) {
    # Propose all the base parameters from the logitnormal rw kinda thing
    # evaluate the proposal probability asap
    par_prop <- c(
      aa = rlogitnorm(n = 1, mu = logit(res_array[ii - 1, 1, "aa"]), sigma = prop_sigma["aa"]),
      bb = rlogitnorm(n = 1, mu = logit(res_array[ii - 1, 1, "bb"]), sigma = prop_sigma["bb"]),
      cc = rlogitnorm(n = 1, mu = logit(res_array[ii - 1, 1, "cc"]), sigma = prop_sigma["cc"]),
      dd = rlogitnorm(n = 1, mu = logit(res_array[ii - 1, 1, "dd"]), sigma = prop_sigma["dd"]),
      ee = rlogitnorm(n = 1, mu = logit(res_array[ii - 1, 1, "ee"]), sigma = prop_sigma["ee"]),
      ff = rlogitnorm(n = 1, mu = logit(res_array[ii - 1, 1, "ff"]), sigma = prop_sigma["ff"]),
      gg = rlogitnorm(n = 1, mu = logit(res_array[ii - 1, 1, "gg"]), sigma = prop_sigma["gg"]),
      hh = rlogitnorm(n = 1, mu = logit(res_array[ii - 1, 1, "hh"]), sigma = prop_sigma["hh"]),
      ww = rlogitnorm(n = 1, mu = logit(res_array[ii - 1, 1, "ww"]), sigma = prop_sigma["ww"])
    )

    # log(q(current | prop))
    log_q_curr_prop <- sum(dlogitnorm(
      x = res_array[ii - 1, 1, ],
      mu = logit(par_prop),
      sigma = prop_sigma,
      log = TRUE
    ))

    # log(q(prop | current))
    log_q_prop_curr <- sum(dlogitnorm(
      x = par_prop,
      mu = logit(res_array[ii - 1, 1, ]),
      sigma = prop_sigma,
      log = TRUE
    ))    

    # use the log_prob method of the stan object (translate to upars first)
    ## Make sure we use the right stan model
    log_probability_prop <- rstan::log_prob(
      model_fit, 
      upars = unconstrain_pars(
        model_fit,
        pars = as.list(par_prop)
      )
    )

    log_probability_curr <- rstan::log_prob(
      model_fit,
      upars = unconstrain_pars(
        model_fit,
        pars = as.list(res_array[ii - 1, 1, ])
      )
    ) 

    # use the wsre est for the prior marginal
    # compute p_12 prop and p_12 current
    with(as.list(par_prop), {
      db <- dd * bb
      e1ab <- ee * (1 - aa - bb)
      p12_prop <<- (db + ww * e1ab) / (db + e1ab)
    })

    with(as.list(res_array[ii - 1, 1, ]), {
      db <- dd * bb
      e1ab <- ee * (1 - aa - bb)
      p12_curr <<- (db + ww * e1ab) / (db + e1ab)
    })

    log_wsre_term <- log(evaluate(
      wsre_obj = wsre_est,
      x_nu = p12_curr,
      x_de = p12_prop 
    ))

    # usual MH correction - 
    log_alpha <- 
      log_probability_prop - log_probability_curr +
      log_q_curr_prop - log_q_prop_curr +
      log_wsre_term

    if (runif(1) < exp(log_alpha)) {
      res_array[ii, 1, ] <- par_prop 
    } else {
      res_array[ii, 1, ] <- res_array[ii - 1, 1, ]
    }

    # log for progress - test
    if (ii %% 1000 == 0) {
      flog.info(
        sprintf("Chain: %d, Iteration: %d", chain_id, ii),
        name = "wsre-logger"
      )
    }
  }
  return(res_array)
})

mcmc_samples <- abind::abind(res_list, along = 2)
mcmc_samples <- mcmc_samples[thin_vec, , ]
p1 <- bayesplot::mcmc_trace(mcmc_samples)


saveRDS(
  object = mcmc_samples,
  file = "rds/hiv-example/stage-one-wsre-samples.rds"
)

ggsave_fullwidth(
  "plots/hiv-example/chain-check.pdf",
  p1
)