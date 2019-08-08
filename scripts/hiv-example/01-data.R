library(tibble)
library(dplyr)

hiv_data <- tibble(
  parameter = sprintf("p_%02d", 1:12),
  y = c(11044, 12, 252, 10, 74, 254, 43, 4, 87, 12, 14, 5),
  n = c(104577, 882, 15428, 473, 136139, 102287, 60, 17, 254, 15, 118, 32)
)

saveRDS(
  object = hiv_data,
  file = "rds/hiv-example/hiv-data.rds"
)

pars_of_interest <- c(
  "aa",
  #"zz",
  "bb",
  "cc",
  "dd",
  "ee",
  "ff",
  "gg",
  "hh",
  "ww",
  sprintf("p[%d]", 1:12)
)

saveRDS(
  object = pars_of_interest,
  file = "rds/hiv-example/pars-of-interest.rds"
)