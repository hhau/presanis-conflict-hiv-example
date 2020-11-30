evaluate_telescope_1d <- function(
  wsre_obj,
  x_nu,
  x_de, 
  mc_cores = 1, 
  log_result = FALSE
) {
  # only for '1d' wsre
  stopifnot(length(x_nu) == 1, length(x_de) == 1)
  
  sub_ests <- wsre_obj$estimates
  n_sub_ests <- length(sub_ests)
  
  # step 1 
  # work with x_nu_eval and x_de_eval so we only have to think about this problem
  # in one direction. it is symmetric after all?
  if (x_nu < x_de) {
    x_nu_eval <- x_nu
    x_de_eval <- x_de
    invert_result <- FALSE
  } else {
    x_nu_eval <- x_de
    x_de_eval <- x_nu
    invert_result <- TRUE
  }
  
  # step two
  # figure out which weighted targets have means inbetween x_nu_eval and x_de_eval
  # this could also be a vector, => function is now 1d, think about
  # 2d case later
  emp_means <- sapply(sub_ests, function(x) x$properties$sample_mean)
  emp_means_index <- which((emp_means < x_de_eval) & (x_nu_eval < emp_means))
  n_wfs_to_include <- length(emp_means_index)

  # this is a little global variable mess
  # n_tele_terms_vec[ii - 1] <<- n_wfs_to_include

  # need to check if all emp_means are false (there is no sample_mean 
  # inbetween x_nu_eval and x_de_eval)
  # in which case, find the 'nearest' wf and compute ratio using that?
  if (n_wfs_to_include == 0) {
    point <- mean(c(x_nu_eval, x_de_eval))
    min_l2_dist_wf_index <- which.min((emp_means - point)^2)
    log_ratio <- with(sub_ests[[min_l2_dist_wf_index]], {
      log(ratio(x_nu_eval, x_de_eval))
    })
    if (invert_result) {
      log_ratio <- -1 * log_ratio
    }
    ifelse(log_result, return(log_ratio), return(exp(log_ratio)))
  }
  
  # now we are assuming that length(emp_means_index) > 0
  points_and_index <- lapply(emp_means_index, function(sub_index) {
    with(
      sub_ests[[sub_index]], 
      data.frame(
        point = c(
          properties$sample_lower_quantile, 
          properties$sample_upper_quantile
        ), # need to think about 2d version
        type = c("lower", "upper"),
        wf_index = sub_index
      )
    )    
  })
  
  # if there is only 1 point, compute ratio using appropriate wf and move on
  if (n_wfs_to_include == 1) {
    index <- points_and_index[[1]]$wf_index[1]
    log_ratio <- with(sub_ests[[index]], {
      log(ratio(x_nu_eval, x_de_eval))
    })
    if (invert_result) {
      log_ratio <- -1 * log_ratio
    }
    ifelse(log_result, return(log_ratio), return(exp(log_ratio)))
  }
  
  # # step three
  # # check to see if the smallest 10th quantile is less than x_nu_eval, if it is, remove it
  min_index <- points_and_index[[1]]$wf_index[1]
  max_index <- points_and_index[[n_wfs_to_include]]$wf_index[1]
  
  min_tele_point <- points_and_index[[1]]$point[1]
  max_tele_point <- points_and_index[[n_wfs_to_include]]$point[2]
  
  if (x_nu_eval > min_tele_point) {
    points_and_index[[1]] <- points_and_index[[1]][-1, ]
  }
  
  # # check to see if the largest 10th quantile is greater than x_de_eval, if it is, remove it
  if (x_de_eval < max_tele_point) {
    points_and_index[[n_wfs_to_include]] <- points_and_index[[n_wfs_to_include]][-2, ]
  }
  
  
  # # we should now be able to form a df like this:
  # point | type | wf_index 
  all_point_df <- do.call(rbind, c(
    list(data.frame(point = x_nu_eval, type = "eval", wf_index = min_index)),
    points_and_index,
    list(data.frame(point = x_de_eval, type = "eval", wf_index = max_index)),
    make.row.names = FALSE
  ))
  
  # sprintf("x_nu: %f, x_de: %f", x_nu, x_de)
  # print(all_point_df)
  
  # # these are the telescoping quantities
  log_ratio_list <- list()
  
  # this for-loop could be parallel-lise-able
  for (ii in 2 : nrow(all_point_df)) {
    local_nu <- all_point_df[ii - 1, 1]
    local_de <- all_point_df[ii, 1]
    nu_index <- all_point_df[ii - 1, 3]
    de_index <- all_point_df[ii, 3]
    
    if (nu_index == de_index) {
      # compute using sensible index
      log_ratio_list[[ii - 1]] <- with(sub_ests[[nu_index]], {
        log(ratio(local_nu, local_de))
      })
    } else {
      nu_est <- sub_ests[[nu_index]]
      weight_nu <- nu_est$weighting(local_nu, local_de)
      ratio_nu <- nu_est$ratio(local_nu, local_de)
      prod_nu <- weight_nu * ratio_nu

      if (!wsre:::.is_numerically_okay(prod_nu)) {
        stop("prod_nu not numerically okay")
      }
      
      de_est <- sub_ests[[de_index]]
      weight_de <- de_est$weighting(local_nu, local_de)
      ratio_de <- de_est$ratio(local_nu, local_de)
      prod_de <- weight_de * ratio_de

      if (!wsre:::.is_numerically_okay(prod_de)) {
        stop("prod_de not numerically okay")
      }

      log_ratio_list[[ii - 1]] <- -log(weight_nu + weight_de) + 
        log(prod_nu + prod_de)
    } 
  }

  log_res_sum <- sum(unlist(log_ratio_list))
  if (invert_result) {
    log_res_sum <- -1 * log_res_sum
  }
  
  ifelse(log_result, return(log_res_sum), return(exp(log_res_sum)))
}

evaluate_telescope_1d_means <- function(
  wsre_obj,
  x_nu,
  x_de, 
  mc_cores = 1, 
  log_result = FALSE
) {
  # only for '1d' wsre
  stopifnot(length(x_nu) == 1, length(x_de) == 1)
  
  sub_ests <- wsre_obj$estimates
  n_sub_ests <- length(sub_ests)
  
  # step 1 
  # work with x_nu_eval and x_de_eval so we only have to think about this problem
  # in one direction. it is symmetric after all?
  if (x_nu < x_de) {
    x_nu_eval <- x_nu
    x_de_eval <- x_de
    invert_result <- FALSE
  } else {
    x_nu_eval <- x_de
    x_de_eval <- x_nu
    invert_result <- TRUE
  }
  
  # step two
  # figure out which weighted targets have means inbetween x_nu_eval and x_de_eval
  # this could also be a vector, => function is now 1d, think about
  # 2d case later
  emp_means <- sapply(sub_ests, function(x) x$properties$sample_mean)
  emp_means_index <- which((emp_means < x_de_eval) & (x_nu_eval < emp_means))
  n_wfs_to_include <- length(emp_means_index)

  # this is a little global variable mess
  # n_tele_terms_vec[ii - 1] <<- n_wfs_to_include

  # need to check if all emp_means are false (there is no sample_mean 
  # inbetween x_nu_eval and x_de_eval)
  # in which case, find the 'nearest' wf and compute ratio using that?
  if (n_wfs_to_include == 0) {
    point <- mean(c(x_nu_eval, x_de_eval))
    min_l2_dist_wf_index <- which.min((emp_means - point)^2)
    log_ratio <- with(sub_ests[[min_l2_dist_wf_index]], {
      log(ratio(x_nu_eval, x_de_eval))
    })
    if (invert_result) {
      log_ratio <- -1 * log_ratio
    }
    ifelse(log_result, return(log_ratio), return(exp(log_ratio)))
  }
  
  # now we are assuming that length(emp_means_index) > 0
  points_and_index <- lapply(emp_means_index, function(sub_index) {
    with(
      sub_ests[[sub_index]], 
      data.frame(
        point = c(
          properties$sample_mean
        ), # need to think about 2d version
        type = c("mean"),
        wf_index = sub_index
      )
    )    
  })
  
  # if there is only 1 point, compute ratio using appropriate wf and move on
  if (n_wfs_to_include == 1) {
    index <- points_and_index[[1]]$wf_index[1]
    log_ratio <- with(sub_ests[[index]], {
      log(ratio(x_nu_eval, x_de_eval))
    })
    if (invert_result) {
      log_ratio <- -1 * log_ratio
    }
    ifelse(log_result, return(log_ratio), return(exp(log_ratio)))
  }
  
  # # step three
  # # check to see if the smallest 10th quantile is less than x_nu_eval, if it is, remove it
  min_index <- points_and_index[[1]]$wf_index[1]
  max_index <- points_and_index[[n_wfs_to_include]]$wf_index[1]
  
  # min_tele_point <- points_and_index[[1]]$point[1]
  # max_tele_point <- points_and_index[[n_wfs_to_include]]$point[2]
  
  # if (x_nu_eval > min_tele_point) {
  #   points_and_index[[1]] <- points_and_index[[1]][-1, ]
  # }
  
  # # # check to see if the largest 10th quantile is greater than x_de_eval, if it is, remove it
  # if (x_de_eval < max_tele_point) {
  #   points_and_index[[n_wfs_to_include]] <- points_and_index[[n_wfs_to_include]][-2, ]
  # }
  
  
  # # we should now be able to form a df like this:
  # point | type | wf_index 
  all_point_df <- do.call(rbind, c(
    list(data.frame(point = x_nu_eval, type = "eval", wf_index = min_index)),
    points_and_index,
    list(data.frame(point = x_de_eval, type = "eval", wf_index = max_index)),
    make.row.names = FALSE
  ))
  
  # sprintf("x_nu: %f, x_de: %f", x_nu, x_de)
  # print(all_point_df)
  
  # # these are the telescoping quantities
  log_ratio_list <- list()
  
  # this for-loop could be parallel-lise-able
  for (ii in 2 : nrow(all_point_df)) {
    local_nu <- all_point_df[ii - 1, 1]
    local_de <- all_point_df[ii, 1]
    nu_index <- all_point_df[ii - 1, 3]
    de_index <- all_point_df[ii, 3]
    
    if (nu_index == de_index) {
      # compute using sensible index
      log_ratio_list[[ii - 1]] <- with(sub_ests[[nu_index]], {
        log(ratio(local_nu, local_de))
      })
    } else {
      nu_est <- sub_ests[[nu_index]]
      weight_nu <- nu_est$weighting(local_nu, local_de)
      ratio_nu <- nu_est$ratio(local_nu, local_de)
      prod_nu <- weight_nu * ratio_nu

      if (!wsre:::.is_numerically_okay(prod_nu)) {
        stop("prod_nu not numerically okay")
      }
      
      de_est <- sub_ests[[de_index]]
      weight_de <- de_est$weighting(local_nu, local_de)
      ratio_de <- de_est$ratio(local_nu, local_de)
      prod_de <- weight_de * ratio_de

      if (!wsre:::.is_numerically_okay(prod_de)) {
        stop("prod_de not numerically okay")
      }

      log_ratio_list[[ii - 1]] <- -log(weight_nu + weight_de) + 
        log(prod_nu + prod_de)
    } 
  }

  log_res_sum <- sum(unlist(log_ratio_list))
  if (invert_result) {
    log_res_sum <- -1 * log_res_sum
  }
  
  ifelse(log_result, return(log_res_sum), return(exp(log_res_sum)))
}  

# evaluate_telescope_1d_means(wsre_estimate, 0.25, 0.4)

# =============================================================================
# 
# 
evaluate_telescope_1d_dumb <- function(
  wsre_obj,
  x_nu,
  x_de, 
  mc_cores = 1, 
  log_result = FALSE
) {
  # only for '1d' wsre
  stopifnot(length(x_nu) == 1, length(x_de) == 1)
  
  sub_ests <- wsre_obj$estimates
  n_sub_ests <- length(sub_ests)
  
  # step 1 
  # work with x_nu_eval and x_de_eval so we only have to think about this problem
  # in one direction. it is symmetric after all?
  if (x_nu < x_de) {
    x_nu_eval <- x_nu
    x_de_eval <- x_de
    invert_result <- FALSE
  } else {
    x_nu_eval <- x_de
    x_de_eval <- x_nu
    invert_result <- TRUE
  }
  
  # step two
  # figure out which weighted targets have means inbetween x_nu_eval and x_de_eval
  # this could also be a vector, => function is now 1d, think about
  # 2d case later
  emp_means <- sapply(sub_ests, function(x) x$properties$sample_mean)
  emp_means_index <- which((emp_means < x_de_eval) & (x_nu_eval < emp_means))
  n_wfs_to_include <- length(emp_means_index)

  # this is a little global variable mess
  # n_tele_terms_vec[ii - 1] <<- n_wfs_to_include

  # need to check if all emp_means are false (there is no sample_mean 
  # inbetween x_nu_eval and x_de_eval)
  # in which case, find the 'nearest' wf and compute ratio using that?
  if (n_wfs_to_include == 0 || n_wfs_to_include == 1) {
    log_ratio <- log(wsre::evaluate(wsre_obj, x_nu_eval, x_de_eval))
    if (invert_result) {
      log_ratio <- -1 * log_ratio
    }
    ifelse(log_result, return(log_ratio), return(exp(log_ratio)))
  }
  
  # now we are assuming that length(emp_means_index) > 0
  points_and_index <- lapply(emp_means_index, function(sub_index) {
    with(
      sub_ests[[sub_index]], 
      data.frame(
        point = c(
          properties$sample_lower_quantile, 
          properties$sample_upper_quantile
        ), # need to think about 2d version
        type = c("lower", "upper"),
        wf_index = sub_index
      )
    )    
  })
  
  # # step three
  # # check to see if the smallest 10th quantile is less than x_nu_eval, if it is, remove it
  min_index <- points_and_index[[1]]$wf_index[1]
  max_index <- points_and_index[[n_wfs_to_include]]$wf_index[1]
  
  min_tele_point <- points_and_index[[1]]$point[1]
  max_tele_point <- points_and_index[[n_wfs_to_include]]$point[2]
  
  if (x_nu_eval > min_tele_point) {
    points_and_index[[1]] <- points_and_index[[1]][-1, ]
  }
  
  # # check to see if the largest 10th quantile is greater than x_de_eval, if it is, remove it
  if (x_de_eval < max_tele_point) {
    points_and_index[[n_wfs_to_include]] <- points_and_index[[n_wfs_to_include]][-2, ]
  }

  # # we should now be able to form a df like this:
  # point | type | wf_index 
  all_point_df <- do.call(rbind, c(
    list(data.frame(point = x_nu_eval, type = "eval", wf_index = min_index)),
    points_and_index,
    list(data.frame(point = x_de_eval, type = "eval", wf_index = max_index)),
    make.row.names = FALSE
  ))
  
  # # these are the telescoping quantities
  log_ratio_list <- list()
  
  # this for-loop could be parallel-lise-able
  for (ii in 2 : nrow(all_point_df)) {
    local_nu <- all_point_df[ii - 1, 1]
    local_de <- all_point_df[ii, 1]
    nu_index <- all_point_df[ii - 1, 3]
    de_index <- all_point_df[ii, 3]

    log_ratio_list[[ii - 1]] <- log(wsre::evaluate(
      wsre_obj,
      local_nu,
      local_de
    ))
  }

  log_res_sum <- sum(unlist(log_ratio_list))
  if (invert_result) {
    log_res_sum <- -1 * log_res_sum
  }
  
  ifelse(log_result, return(log_res_sum), return(exp(log_res_sum)))
}

evaluate_telescope_fixed_N <- function(
  wsre_obj,
  x_nu,
  x_de,
  mc_cores = 1,
  N_tele_terms = 5,
  log_result = FALSE
) {
  # this should also work for 2D
  sub_ests <- wsre_obj$estimates
  n_sub_ests <- length(sub_ests)
  
  stopifnot(
    length(x_nu) == length(sub_ests[[1]]$properties$sample_mean),
    length(x_de) == length(x_nu)
  )
  
  # how many dimensions
  n_dimensions <- length(x_nu)
  
  # N_tele_terms length list of n_dim arrays
  matrix_of_points <- lapply(1 : n_dimensions, function(d) {
    seq(from = x_nu[d], to = x_de[d], length.out = N_tele_terms)
  }) %>% do.call(cbind, .)
    
  # evaluate all the terms
  log_wsre_terms <- list()
  
  for (ii in 2 : N_tele_terms) {
    log_wsre_terms[[ii - 1]] <- log(wsre::evaluate(
      wsre_obj,
      x_nu = matrix_of_points[ii - 1, ],
      x_de = matrix_of_points[ii, ]
    ))
  }
  
  res <- sum(unlist(log_wsre_terms))
  ifelse(log_result, return(res), return(exp(res)))
}

evaluate_telescope_fixed_dist <- function(
  wsre_obj,
  x_nu,
  x_de,
  mc_cores = 1,
  tele_dist = numeric(), # length(x_nu) vector of distances
  log_result = FALSE
) {
  # this should also work for 2D
  sub_ests <- wsre_obj$estimates
  n_sub_ests <- length(sub_ests)
  
  stopifnot(
    length(x_nu) == length(sub_ests[[1]]$properties$sample_mean),
    length(x_de) == length(x_nu)
  )
  
  if (!missing(tele_dist)) {
    stopifnot(length(tele_dist) == length(x_de))
  }
  
  # how many dimensions
  n_dimensions <- length(x_nu)
  
  # how far to move
  tele_dist <- lapply(sub_ests, function(x) x$properties$sample_sd) %>% 
    do.call(rbind, .) %>% 
    apply(2, mean)
  
  list_of_possible_points <- lapply(1 : n_dimensions, function(d) {
    local_by <- ifelse(x_nu[d] < x_de[d], tele_dist[d], -tele_dist[d])
    temp_seq <- seq(from = x_nu[d], to = x_de[d], by = local_by)
    c(temp_seq, x_de[d])
  }) 
  
  # find the list element with the maximum number of points, replicate
  # the last element of the shorter list to match the number of elements in 
  # the longer list
  max_sub_list_length <- lapply(list_of_possible_points, function(x) length(x)) %>% 
    unlist() %>% 
    max()
  
  matrix_of_points <- lapply(list_of_possible_points, function(x) {
    if (length(x) == max_sub_list_length) {
      return(x)
    }
    n <- length(x)
    res <- c(x, rep(x[n], max_sub_list_length - n))
  }) %>% 
    do.call(cbind, .)
  
  ## TODO: figure out what it means to move a fixed distance when the 
  ## bases differ in scales?
  
  log_wsre_terms <- list()
  
  for (ii in 2 : nrow(matrix_of_points)) {
    log_wsre_terms[[ii - 1]] <- log(wsre::evaluate(
      wsre_obj,
      x_nu = matrix_of_points[ii - 1, ],
      x_de = matrix_of_points[ii, ]
    ))
  }
  
  res <- sum(unlist(log_wsre_terms))
  ifelse(log_result, return(res), return(exp(res)))
}


# new test
# x_de <- 0.8
# x_nu <- 0.2
# evaluate_telescope_1d(wsre_estimate, x_nu, x_de, log_result = FALSE)
# wsre::evaluate(wsre_estimate, x_nu, x_de)
# evaluate_telescope_1d_dumb(wsre_estimate, x_nu, x_de, log_result = FALSE)
# ref_ratio(x_nu, x_de)
# # # random tests - sim random x_nu and x_de and bench / make sure they all work
# # # 
# 
# microbenchmark::microbenchmark(
#   evaluate_telescope_1d(wsre_estimate, x_nu, x_de, log_result = FALSE),
#   wsre::evaluate(wsre_estimate, x_nu, x_de),
#   evaluate_telescope_1d_dumb(wsre_estimate, x_nu, x_de, log_result = FALSE)
#)

# what's going on with the x_nu > 0.75 values?
