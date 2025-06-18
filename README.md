
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MSCCT

<!-- badges: start -->
<!-- badges: end -->

This package contains tests for comparison or two or more survival
curves when the proportional hazards hypothesis is not verified, in
particular when the survival curves cross each other.

## Installation

You can install the development version of MSCCT from
[GitHub](https://github.com/) with:

``` r
devtools::install_github("https://github.com/HMinP/MSCCT")
```

## Example

``` r
library(MSCCT)
```

This package contains:

- The weighted log-rank test

The log-rank test compares for each group and for each time of event the
expected and the observed number of events. The weighted log-rank adds
weights to each time of event. Some (implemented) exemples are the
Flemming-Harrington test and the Gehan-Wilcoxon test. It is also
possible to chose the weights you want.

``` r
multiLR(data_under_PH)
#> (Multiple) Weighted log-rank test 
#> 
#> Weighting : Classic log-rank test 
#> Degrees of freedom : 2 
#> 
#>        Statistic p
#> Test 1  80.17764 0
```

``` r
multiLR(data_under_PH, weights=numeric(), test="fh", rho=1, gamma=0)
#> (Multiple) Weighted log-rank test 
#> 
#> Weighting : Flemming-Harrington test 
#> Parameters : rho = 1 , gamma =  0 
#> Degrees of freedom : 2 
#> 
#>        Statistic            p
#> Test 1  71.47735 3.330669e-16
```

- The Restricted Mean Survival Test

The Restricted Mean Survival Time at time $\tau$ is the area under a
survival curve up to time $\tau$. The RMST test compares the areas under
the survival curves and tests the equality to zero of the difference of
RMST.

``` r
multirmst(data_under_PH, tau=12)
#> (Multiple) test of RMST 
#> Truncation time : 12  
#> Correction : bonferroni 
#> 
#>            dRMST        sd            p   p adjusted
#> 0 VS 1 -1.518564 0.3919589 1.069344e-04 3.208033e-04
#> 0 VS 2 -2.621846 0.4185860 3.762464e-10 1.128739e-09
#> 1 VS 2 -1.103282 0.4461470 1.340176e-02 4.020527e-02
#>  
#> p=1.128739e-09
```
