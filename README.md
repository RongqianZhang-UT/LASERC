
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

The following shows a scheme which shows you how to use `LASERC()` for ABIDE data. Since we use warm start for sequential initialization, we have three modes to run `LASERC()`: (i) eigen initialization for the unpenalized EM phase, (ii) warm-start initialization for the unpenalized EM phase, and (iii) warm-start initialization for the TLP phase. For example, we start with `L=2` with eigen initialization, and use estimates of `L=2` from unpenalized EM phase to initialize unpenalized EM phase of `L=6` for a warm start and use estimates of `L=6` from unpenalized EM phase to initialize TLP phase of `L=6`. Repeat this pattern as `L` increases.

``` r
library(LASERC)

load('data.RData')


### L=2 with eigen initialization in the unpenalized phase

results_L2_OLS <- LASERC(
  dat = Yraw,
  batch = grp,
  mod = x[,1:7],
  L = 2,
  penalty='OLS',initial='eigen',Theta_ini = NULL
  
)

### L=6 with warm start initialization in the unpenalized phase

##initial estimates from L=2 of unpenalized phase

initial_warmL=list(A=results_L2_OLS$estimates$A,S=results_L2_OLS$estimates$S,U=results_L2_OLS$estimates$U,Fmat=results_L2_OLS$estimates$Fmat,xFmat=results_L2_OLS$estimates$xFmat,sigma2=results_L2_OLS$estimates$sigma2)

results_L6_OLS <- LASERC(
  dat = Yraw,
  batch = grp,
  mod = x[,1:7],
  L = 6,
  penalty='OLS+warmL',initial='warmL',Theta_ini = initial_warmL,
  MaxIteration = 100
  
)

### L=6 in the TLP phase

initial_OLS=list(A=results_L6_OLS$estimates$A,S=results_L6_OLS$estimates$S,U=results_L6_OLS$estimates$U,xFmat=results_L6_OLS$estimates$xFmat,Fmat=results_L6_OLS$estimates$Fmat,sigma2=results_L6_OLS$estimates$sigma2,phi2=results_L6_OLS$estimates$phi2)

results_LASERC <- LASERC(
  dat = Yraw,
  batch = grp,
  mod = x[,1:7],
  L = 6,
  penalty='TLP',initial='warm',Theta_ini = initial_OLS,
  MaxIteration = 100
  
)

```

