
# LASERC

<!-- badges: start -->
<!-- badges: end -->

The goal of LASERC is to provides functions for site-effect correction and data harmonization for functional connectivity data in multi-site fMRI studies using a latent factor model. 

## Installation

You can install the development version of LASERC like so:

``` r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("RongqianZhang-UT/LASERC")
```

## Example

This is a basic example which shows you how to use `LASERC()`:

``` r
library(LASERC)

data(sim_data)
dat=sim_data$dat
batch=sim_data$batch
mod=sim_data$mod
n=nrow(dat)
p=ncol(dat)
V=(sqrt(1+8*p)-1)/2
results <- LASERC(
  dat = dat,
  batch = batch,
  mod = mod,
  L = 5
)
```

