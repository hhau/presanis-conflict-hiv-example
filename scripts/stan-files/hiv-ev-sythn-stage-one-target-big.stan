functions {
  real log_kde(real x, vector x_samples, real bandwidth) {
    int n_samples = num_elements(x_samples);
    vector [n_samples] temp_log_vector;
    real res;

    for (ii in 1 : n_samples) {
      temp_log_vector[ii] = normal_lpdf((x - x_samples[ii]) / bandwidth | 0.0, 1.0);
    } 

    res = log_sum_exp(temp_log_vector);
    return(res);
  }
}

data {
  int n_studies;
  int y_obs [n_studies];
  int n_obs [n_studies];

  // KDE parameters
  int n_prior_samples;
  vector [n_prior_samples] phi_prior_samples;
  real <lower = 0> bandwidth;
}

parameters {
  real <lower = 0, upper = 1> aa;
  real <lower = 0, upper = 1> bb;
  real <lower = 0, upper = 1> cc;
  real <lower = 0, upper = 1> dd;
  real <lower = 0, upper = 1> ee;
  real <lower = 0, upper = 1> ff;
  real <lower = 0, upper = 1> gg;
  real <lower = 0, upper = 1> hh;
  real <lower = 0, upper = 1> ww;
}

transformed parameters {
  real <lower = 0, upper = 1> p [n_studies];

  {
    // doubles
    real ca = cc * aa;
    real db = dd * bb;

    // triples
    real e1ab = ee * (1 - aa - bb);
    real fca = ff * cc * aa;
    real gdb = gg * dd * bb;

    // quadruples
    real he1ab = hh * e1ab;
    real we1ab = ww * e1ab;

    // probabilities
    p[1] = aa;
    p[2] = bb;
    p[3] = cc;
    p[4] = dd;
    p[5] = (db + e1ab) / (1 - aa);
    p[6] = ca + db + e1ab;
    p[7] = (fca) / (fca + gdb + he1ab);
    p[8] = (gdb) / (gdb + he1ab);
    p[9] = (fca + gdb + he1ab) / (ca + db + e1ab);
    p[10] = gg;
    p[11] = ww;
    p[12] = (db + we1ab) / (db + e1ab);  
  }
}

model {
  // jeffereys priors
  aa ~ beta(1, 2);
  bb ~ beta(1, 9);
  cc ~ beta(1, 9);
  dd ~ beta(1, 9);
  ee ~ beta(1, 9);
  ff ~ beta(1, 1);
  gg ~ beta(1, 1);
  hh ~ beta(1, 1);
  ww ~ beta(3, 1);

  // likelihood - drops the 12th study
  y_obs[1 : (n_studies - 1)] ~ binomial(n_obs[1 : (n_studies - 1)], p[1 : (n_studies - 1)]);

  // marginalise out the prior on phi
  target += -1 * log_kde(p[12], phi_prior_samples, bandwidth);
}
