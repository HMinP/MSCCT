
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MSCCT

<!-- badges: start -->
<!-- badges: end -->

Multiple Survival Crossing Curves Tests

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
multi_lr(data_under_PH)
#> (Multiple) Weighted log-rank test 
#> 
#> Weighting : Classic log-rank test 
#> Degrees of freedom : 2 
#> 
#>        Statistic p
#> Test 1  80.17764 0
```

``` r
multi_lr(data_under_PH, weights=numeric(), test="fh", rho=1, gamma=0)
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
multi_rmst(data_under_PH, tau=12, nboot=100, method="bonferroni")
#> (Multiple) test of RMST 
#> Truncation time : 12  
#> Correction : bonferroni 
#> 
#>            dRMST        sd            p   p adjusted
#> 0 VS 1 -1.518564 0.4162182 2.637962e-04 7.913886e-04
#> 0 VS 2 -2.621846 0.4092926 1.495843e-10 4.487530e-10
#> 1 VS 2 -1.103282 0.4060768 6.589085e-03 1.976726e-02
#>  
#> p=4.48753e-10
```

- The Two-stage test

The two-stage test is a combination of two tests. The first one is a
classic log-rank test. When the log-rank test is not significant, this
means that the survival curves are either different or they cross each
other and the log-rank is not powerful enough. In order to differentiate
these cases, a second test is performed. This second test is a weighted
log-rank test with weights that allows to differentiate the two previous
cases.

``` r
multi_ts(data_under_PH, eps=0.1, nboot=100, method="BH")
#> (Multiple) Two-Staged test 
#> Correction : BH  
#> 
#>                  p1   p2            p        adj_p
#> 0 VS 1 9.377383e-09 0.75 1.390618e-07 2.085926e-07
#> 0 VS 2 0.000000e+00 0.00 3.774758e-15 1.132427e-14
#> 1 VS 2 6.246677e-04 0.14 9.046541e-04 9.046541e-04
#>  
#> p=1.132427e-14
```
