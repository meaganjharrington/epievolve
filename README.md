
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rodent package

<!-- badges: start -->

<!-- badges: end -->

The goal of rodent is to estimate time-varying Rt(t) using the
odin-monty framework. Estimates Rt(t) from time series of case counts
using a mechanistic SIR model. Enables interpretable inference on
epidemic dynamics and intervention effects. \## Installation You can
install the development version of rodent from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("meaganjharrington/rodent")
```

## Example: 2009 Influenza A(H1N1), Pennsylvania school outbreak

This example uses the `Flu2009` dataset from the `EpiEstim` package. The
outbreak comprised 129 cases over 32 days in a school of 369 students,
with a known school closure from day 18 to 24.

Change points are set at days 1, 8, 18, and 25, with the third block
capturing the period of school closure.

``` r
library(rodent)
library(EpiEstim)
library(dplyr)

# Prepare incidence data
data(Flu2009)
flu_daily <- Flu2009$incidence |>
  mutate(time = row_number()) |>
  rename(cases = I) |>
  select(time, cases) |>
  as.data.frame()

# Fit rodent
out_flu <- estimate_R_blocks(
  data        = flu_daily,
  N           = 369,
  gamma       = 1 / 2.6,
  beta_breaks = c(1, 8, 18, 25),
  priors = list(mean_beta = log(0.5), sd_beta = 0.5,
                rw_sd_beta = 0.3, mean_I0 = log(5), sd_I0 = 0.5),
  inits  = list(beta = 0.5, I0 = 5),
  mcmc   = list(n_steps = 15000, burnin = 0.5,
                seed = 444, n_chains = 3)
)

# View Rt series and plot
head(out_flu$Rt_series)
out_flu$Rt_plot
```

Blockwise Rt estimates decline over the outbreak, from ~1.36 in the
first block to ~0.86 in the final block, with a notable drop to ~0.41
during the school closure period (days 18–24).
