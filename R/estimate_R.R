#' Package Workflow Overview
#'
#' High-level idea:
#'   1) Model transmission with an SIR system where β(t) is piecewise-constant
#'      over user-defined time blocks, γ is fixed, and the initial infectious count I0
#'      is unknown. R_t is computed as R_t = (β(t) / γ) * S(t) / N
#'   2) Infer the unknown parameters (β-blocks and I0) by sampling a posterior
#'      using a random-walk MH sampler. The likelihood is evaluated by integrating
#'      the SIR system against the observed daily incidence
#'   3) Propagate posterior uncertainty by simulating S(t) for many posterior
#'      draws and computing pointwise R_t quantiles
#'
#' Model specification
#' - States: S(t), I(t), R(t) with S + I + R = N.
#' - Transmission: β(t) = β_k for t in block k (piecewise constant).
#' - Removal/recovery: fixed γ.
#' - Initial conditions: I0 unknown (inferred), S0 = N - I0 (no pre-existing immunity
#'   in the basic version; can be extended to include R0 > 0 if needed).
#' - Reproduction number:
#'     R_t = (β(t) / γ) * S(t) / N
#'
#' Rt(t) with stepwise beta(t), fixed gamma; infer I0 and beta blocks via Metropolis-Hastings random walk MCMC
#' Using full odin2 -> dust2 -> monty workflow (deterministic likelihood)
#' log-posterior = loglikelihood + logprior
#'
#' @param incidence data.frame with columns: time (consecutive integers), cases >= 0.
#' Time can start at any integer, we keep it and set time_start to "first time - 1".
#'   - time : consecutive integers (any starting value allowed)
#'   - cases: non-negative integers (daily incidence)
#' @param N Numeric scalar > 0, population size
#' @param gamma Numeric scalar > 0, fixed recovery rate
#' R Immunity existing in population, default = 0 - add back in once baseline immunity functionality added!
#' @param beta_breaks Integer vector in ORIGINAL time units; mapped to indices internally
#' @param mcmc List with settings:
#'                - n_steps (6000)
#'                - burnin (0.5)
#'                - proposal (NULL -> diag(0.02^2))
#'                - seed (4)
#'                - n_rt_draws (300) posterior draws for Rt uncertainty bands
#' @param priors List: mean_beta=log(0.3), sd_beta=0.5, rw_sd_beta=0.2,
#'               mean_I0=log(10), sd_I0=0.5
#' @param inits  List: beta=0.3, I0=10
#' @return list(samples, estimates, Rt_series, fit, diagnostics, model_used, blocks, fixed)
#' @export
estimate_R_blocks <- function(data,
                              N,
                              gamma,
                              beta_breaks,
                              priors = list(mean_beta  = log(0.3),
                                            sd_beta    = 0.5,
                                            rw_sd_beta = 0.2,
                                            mean_I0    = log(10),
                                            sd_I0      = 0.5),
                              inits = list(beta = 0.3, I0 = 10),
                              mcmc  = list(n_steps  = 6000,
                                           n_chains = 3,
                                           burnin   = 0.5,
                                           seed     = 4)) {

  ## 0. Setup
  data  <- data[order(data$time), , drop = FALSE]
  time  <- data$time
  Tlen  <- length(time)

  time_start <- min(time) - 1L

  # beta grid to cover model start for interpolation to work
  beta_times <- sort(unique(beta_breaks))
  if (beta_times[1] > time_start) {
    beta_times <- c(time_start, beta_times)
  }
  K <- length(beta_times)

  ## 1. odin2 model
  # interpolate() for beta(t)
  sir <- make_sir_beta_blocks_generator() # see other file

  ## 2. dust unfilter
  filter <- dust2::dust_unfilter_create(
    sir,
    data       = data,
    time_start = time_start
  )

  ## 3. Parameter names
  # beta_1..beta_K, I0
  par_names <- c(sprintf("beta_%d", seq_len(K)), "I0")

  ## 4. Likelihood via dust_likelihood_run()
  likelihood <- monty::monty_model(
    list(
      parameters = par_names,
      density = function(theta) {
        beta_vals <- theta[seq_len(K)]
        I0_val    <- theta[K + 1L]

        if (any(beta_vals <= 0) || I0_val <= 0) return(-Inf)

        pars <- list(
          N           = N,
          gamma       = gamma,
          I0          = I0_val,
          n_beta      = K,
          beta_times  = beta_times,
          beta_values = beta_vals
        )

        dust2::dust_likelihood_run(filter, pars)
      }
    )
  )

  ## 5. Prior beta_k and I0
  # random-walk on log beta
  prior <- monty::monty_model(
    list(
      parameters = par_names,
      density = function(theta) {
        beta_vals <- theta[seq_len(K)]
        I0_val    <- theta[K + 1L]

        if (any(beta_vals <= 0) || I0_val <= 0) return(-Inf)

        log_beta <- log(beta_vals)

        lp <- sum(dnorm(log_beta,
                        mean = priors$mean_beta,
                        sd   = priors$sd_beta,
                        log  = TRUE))

        if (K > 1)
          lp <- lp + sum(dnorm(diff(log_beta),
                               mean = 0,
                               sd   = priors$rw_sd_beta,
                               log  = TRUE))

        lp + dnorm(log(I0_val),
                   mean = priors$mean_I0,
                   sd   = priors$sd_I0,
                   log  = TRUE)
      }
    )
  )

  ## 6. Posterior
  posterior <- likelihood + prior

  ## 7. Sample
  set.seed(mcmc$seed %||% 1L)

  init_theta <- c(rep(inits$beta, K), inits$I0)

  vcv     <- diag(K + 1L) * 0.05^2
  sampler <- monty::monty_sampler_adaptive(vcv, initial_vcv_weight = 50)

  samples <- monty::monty_sample(
    posterior,
    sampler,
    n_steps  = mcmc$n_steps  %||% 6000L,
    initial  = init_theta,
    n_chains = mcmc$n_chains %||% 3L
  )

  ## 8. Process burnins + extract draws
  burn     <- floor((mcmc$burnin %||% 0.5) * dim(samples$pars)[2])
  thin     <- monty::monty_samples_thin(samples, burnin = burn)
  pars_arr <- thin$pars
  D        <- dim(pars_arr)[2] * dim(pars_arr)[3]
  pars_mat <- matrix(pars_arr, nrow = K + 1L, ncol = D)

  beta_draws <- pars_mat[seq_len(K), , drop = FALSE]
  I0_draws   <- pars_mat[K + 1L, ]

  ## 9. Expand posterior beta(t) for Rt calculation
  expand_beta_blocks <- function(beta_vec, beta_times, tvec) {
    out <- numeric(length(tvec))
    for (k in seq_along(beta_times)) {
      t_start <- beta_times[k]
      t_end   <- if (k < length(beta_times)) beta_times[k + 1] - 1 else max(tvec)
      out[tvec >= t_start & tvec <= t_end] <- beta_vec[k]
    }
    out
  }

  beta_t_mat <- vapply(
    seq_len(D),
    function(d) expand_beta_blocks(beta_draws[, d], beta_times, time),
    numeric(Tlen)
  )

  ## 10. Extract S(t)
  time_start <- min(time) - 1L

  S_mat <- matrix(NA_real_, nrow = Tlen, ncol = D)

  for (d in seq_len(D)) {
    pars_d <- list(
      N           = N,
      gamma       = gamma,
      I0          = I0_draws[d],
      n_beta      = K,
      beta_times  = beta_times,
      beta_values = beta_draws[, d]
    )

    sys <- dust2::dust_system_create(
      sir,
      pars        = pars_d,
      n_particles = 1,
      time        = time_start
    )
    dust2::dust_system_set_state_initial(sys)

    y_sim      <- dust2::dust_system_simulate(sys, time)
    unpacked   <- dust2::dust_unpack_state(sys, y_sim)
    S_mat[, d] <- drop(unpacked$S)
  }

  ## 11. Rt(t) = (beta(t) / gamma) * S(t) / N
  Rt_mat <- (beta_t_mat / gamma) * (S_mat / N)

  qfun <- function(x) stats::quantile(x, c(0.025, 0.5, 0.975))
  fmt  <- function(mat) {
    q <- t(apply(mat, 1, qfun))
    colnames(q) <- c("lwr", "median", "upr")
    data.frame(time = time, q)
  }

  ## 12. Return
  list(
    samples          = thin,
    beta_draws       = beta_draws,
    I0_draws         = I0_draws,
    beta_t_quantiles = fmt(beta_t_mat),
    S_quantiles      = fmt(S_mat),
    Rt_quantiles     = fmt(Rt_mat),
    Rt_mat           = Rt_mat,
    S_mat            = S_mat,
    beta_t_mat       = beta_t_mat,
    beta_times       = beta_times,
    model_components = list(
      gen       = sir,
      filter    = filter,
      posterior = posterior
    )
  )
}
