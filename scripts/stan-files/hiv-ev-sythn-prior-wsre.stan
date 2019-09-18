data {
  int <lower = 1> target_dimension;
  real wf_mean [target_dimension];
  real <lower = 0> wf_sd [target_dimension];
  real <lower = 0> wf_exponent;
   // for future use
}

transformed data {
  int n_studies = 12;  
}

parameters {
  real <lower = 0, upper = 1> aa;
  // real <lower = 0, upper = 1> zz;
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
  // real <lower = 0, upper = 1> bb = zz * (1 - aa);
  real ca = cc * aa;
  real db = dd * bb;
  
  // triples
  real e1ab = ee * (1 - aa - bb);
  real fca = ff * cc * aa;
  real gdb = gg * dd * bb;
  
  // quadruples
  real he1ab = hh * e1ab;
  real we1ab = ww * e1ab;

  // wsre
  real x [target_dimension];
 
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

  x[1] = p[12];
}

model {
  aa ~ beta(1, 2);
  // zz ~ beta(1, 1);
  bb ~ beta(1, 9);
  cc ~ beta(1, 9);
  dd ~ beta(1, 9);
  ee ~ beta(1, 9);
  ff ~ beta(1, 1);
  gg ~ beta(1, 1);
  hh ~ beta(1, 1);
  ww ~ beta(3, 1);
  target += wf_exponent * normal_lpdf(x | wf_mean, wf_sd);
}

generated quantities {
    
}
