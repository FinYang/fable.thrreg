
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fable.thrreg

<!-- badges: start -->
<!-- badges: end -->

The goal of fable.thrreg is to extend `fable` package to fit threshold
regression.

## Installation

You can install the development version of fable.thrreg like so:

``` r
remotes::install_github("FinYang/fable.thrreg")
```

## Example

This is a basic example which shows you how to fit a threshold
regression.

``` r
library(fable.thrreg)

library(tidyverse) 
library(fable)
library(tsibble)
```

``` r
# Simulate a time series
threshold_process <- function(lag1){
  if(abs(lag1) > 1) {
    # Above threshold of 1
    # Mean reverting 
    # AR1 process with coeffcient 0.8
    lag1 * 0.8 + rnorm(1)
  } else {
    # Below threshold of 1
    # a unit root process
    lag1 + rnorm(1)
  }
}
set.seed(2222)
time_span <- 2000
y <- numeric(time_span)
for(i in 2:time_span) y[[i]] <- threshold_process(y[[i-1]])

# Convert to tsibble
df <- tsibble(y = y, idx = seq_along(y), index = idx)
```

``` r

# Fit the threshold regression
fit <- df %>% 
  model(thrreg = THRREG(y ~ offset(lag(y)) + lag(y)*ind(abs(lag(y)) >= gamma(1) )))
# Getting estimates
est <- tidy(fit)
est
#> # A tibble: 3 × 6
#>   .model term                                estimate std.e…¹ statis…²   p.value
#>   <chr>  <chr>                                  <dbl>   <dbl>    <dbl>     <dbl>
#> 1 thrreg (Intercept)                         -0.00136  0.0227  -0.0602  9.52e- 1
#> 2 thrreg I(lag(y) * (abs(lag(y)) >= 0.97949… -0.220    0.0141 -15.6     4.18e-52
#> 3 thrreg .gamma_1                             0.979   NA       NA      NA       
#> # … with abbreviated variable names ¹​std.error, ²​statistic

# Estimated threshold value
.gamma_1 <- est$estimate[est$term == ".gamma_1"] 
# Plot
autoplot(df, y) +
  geom_hline(yintercept = c(-1, 1)) +
  geom_hline(yintercept = c(-.gamma_1, .gamma_1), colour = "red")
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />
