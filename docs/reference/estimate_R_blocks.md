# Package Workflow Overview

High-level idea:

1.  Model transmission with an SIR system where β(t) is
    piecewise-constant over user-defined time blocks, γ is fixed, and
    the initial infectious count I0 is unknown. R_t is computed as R_t =
    (β(t) / γ) \* S(t) / N

2.  Infer the unknown parameters (β-blocks and I0) by sampling a
    posterior using a random-walk MH sampler. The likelihood is
    evaluated by integrating the SIR system against the observed daily
    incidence

3.  Propagate posterior uncertainty by simulating S(t) for many
    posterior draws and computing pointwise R_t quantiles

## Usage

``` r
estimate_R_blocks(
  data,
  N,
  gamma,
  beta_breaks,
  priors = list(mean_beta = log(0.3), sd_beta = 0.5, rw_sd_beta = 0.2, mean_I0 = log(10),
    sd_I0 = 0.5),
  inits = list(beta = 0.3, I0 = 10),
  mcmc = list(n_steps = 6000, n_chains = 3, burnin = 0.5, seed = 4)
)
```

## Arguments

- data:

  data.frame with columns: time (consecutive integers), cases \>= 0.
  Time can start at any integer, we keep it and set time_start to "first
  time - 1".

  - time : consecutive integers (any starting value allowed)

  - cases: non-negative integers (daily incidence)

- N:

  Numeric scalar \> 0, population size

- gamma:

  Numeric scalar \> 0, fixed recovery rate R Immunity existing in
  population, default = 0 - add back in once baseline immunity
  functionality added!

- beta_breaks:

  Integer vector in ORIGINAL time units; mapped to indices internally

- priors:

  List: mean_beta=log(0.3), sd_beta=0.5, rw_sd_beta=0.2,
  mean_I0=log(10), sd_I0=0.5

- inits:

  List: beta=0.3, I0=10

- mcmc:

  List with settings: - n_steps (6000) - burnin (0.5) - proposal (NULL
  -\> diag(0.02^2)) - seed (4) - n_rt_draws (300) posterior draws for Rt
  uncertainty bands

## Value

list(samples, estimates, Rt_series, fit, diagnostics, model_used,
blocks, fixed)

## Details

Model specification

- States: S(t), I(t), R(t) with S + I + R = N.

- Transmission: β(t) = β_k for t in block k (piecewise constant).

- Removal/recovery: fixed γ.

- Initial conditions: I0 unknown (inferred), S0 = N - I0 (no
  pre-existing immunity in the basic version; can be extended to include
  R0 \> 0 if needed).

- Reproduction number: R_t = (β(t) / γ) \* S(t) / N

Rt(t) with stepwise beta(t), fixed gamma; infer I0 and beta blocks via
Metropolis-Hastings random walk MCMC Using full odin2 -\> dust2 -\>
monty workflow (deterministic likelihood) log-posterior =
loglikelihood + logprior
