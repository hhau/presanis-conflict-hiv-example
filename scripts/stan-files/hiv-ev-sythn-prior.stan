data {
  int n_studies;
}

parameters {
  
}

model {

}

generated quantities {
  real <lower = 0, upper = 1> p [n_studies];
  
  // some Ades, Cliff 2002 magic
  real <lower = 0, upper = 1> aa = beta_rng(1, 2);
  // real <lower = 0, upper = 1> zz = beta_rng(1, 1);
  // real <lower = 0, upper = 1> bb = zz * (1 - aa);
  real <lower = 0, upper = 1> bb = beta_rng(1, 9);
  real <lower = 0, upper = 1> cc = beta_rng(1, 9);
  real <lower = 0, upper = 1> dd = beta_rng(1, 9);
  real <lower = 0, upper = 1> ee = beta_rng(1, 9);
  real <lower = 0, upper = 1> ff = beta_rng(1, 1);
  real <lower = 0, upper = 1> gg = beta_rng(1, 1);
  real <lower = 0, upper = 1> hh = beta_rng(1, 1);
  real <lower = 0, upper = 1> ww = beta_rng(3, 1);

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
