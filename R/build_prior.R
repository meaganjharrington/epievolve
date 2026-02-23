#' Build deterministic prior object
#' @keywords internal

build_prior <- function(
  K,
  mean_beta,
  sd_beta,
  rw_sd_beta,
  mean_I0,
  sd_I0,
  loglik,
  beta_max = NULL
  ) {

## Priors (log-scale)
#
  logprior <- function(theta) {
    log_b <- theta[seq_len(K)]; log_I <- theta[K + 1L]

    p_b_ind <- sum(stats::dnorm(log_b, mean = mean_beta,
                                sd = sd_beta, log = TRUE))
    p_b_rw  <- if (K > 1)
      sum(stats::dnorm(diff(log_b), mean = 0,
                       sd = rw_sd_beta, log = TRUE)) else 0
    p_I     <- stats::dnorm(log_I, mean = mean_I0,
                            sd = sd_I0, log = TRUE)

    ## Soft upper wall on beta via max_attack_rate
    ## Penalises any beta block above beta_max using a tight half-normal
    ## The penalty is 0 when beta <= beta_max, increasingly negative above it
    p_attack <- if (!is.null(beta_max)) {
      excess <- log_b - log(beta_max)          # positive when beta > beta_max
      sum(ifelse(excess > 0,
                 stats::dnorm(excess, mean = 0, sd = 0.05, log = TRUE) -
                   stats::dnorm(0,    mean = 0, sd = 0.05, log = TRUE),
                 0))                            # no penalty when below ceiling
    } else 0

    p_b_ind + p_b_rw + p_I + p_attack
  }

## Posterior
density_theta <- function(theta)
  loglik(theta) + logprior(theta) # posterior log density used by monty.

## Return objects needed downstream
list(
  logprior = logprior,
  density_theta = density_theta
)
}
