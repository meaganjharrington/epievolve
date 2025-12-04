#' Simulate SIR Model with Single Time-Varying Transmission Rate
#'
#' @param N Total population size
#' @param I0 Initial number of infected individuals
#' @param gamma Recovery rate (1/infectious period)
#' @param beta1 Transmission rate before change point
#' @param beta2 Transmission rate after change point
#' @param change_time Time point when transmission rate changes
#' @param times Vector of times at which to return output (default: 0 to 100 by 1)
#'
#' @return A list containing:
#'   \item{times}{Time points}
#'   \item{S}{Susceptible counts over time}
#'   \item{I}{Infected counts over time}
#'   \item{R}{Recovered counts over time}
#'   \item{beta}{Transmission rate over time}
#'   \item{parameters}{List of parameters used}
#'
#' @examples
#' # Simple intervention at day 50
#' result <- run_sir_step(
#'   N = 1000,
#'   I0 = 10,
#'   gamma = 0.1,
#'   beta1 = 0.4,
#'   beta2 = 0.1,
#'   change_time = 50,
#'   times = seq(0, 100, by = 0.5)
#' )
#'
#' # Plot results
#' plot(result)
#'
#' @export
run_sir_step <- function(N = 1000,
                         I0 = 10,
                         gamma = 0.1,
                         beta1 = 0.4,
                         beta2 = 0.1,
                         change_time = 50,
                         times = seq(0, 100, by = 1)) {

  # Validate inputs
  if (N <= 0) stop("N must be positive")
  if (I0 <= 0 || I0 >= N) stop("I0 must be positive and less than N")
  if (gamma <= 0) stop("gamma must be positive")
  if (beta1 <= 0 || beta2 <= 0) stop("beta1 and beta2 must be positive")
  if (change_time < min(times)) stop("change_time must be within time range")
  if (!all(diff(times) > 0)) stop("times must be strictly increasing")

  # Define the SIR model
  sir_model <- odin2::odin({
    deriv(S) <- -beta * S * I / N
    deriv(I) <- beta * S * I / N - gamma * I
    deriv(R) <- gamma * I

    initial(S) <- N - I0
    initial(I) <- I0
    initial(R) <- 0

    N           <- parameter(1000)
    I0          <- parameter(10)
    beta        <- if (t < change_time) beta1 else beta2
    beta1       <- parameter(0.4)
    beta2       <- parameter(0.1)
    change_time <- parameter(50)
    gamma       <- parameter(0.1)
  })

  # Create the dust system (deterministic ODE solver)
  sys <- dust2::dust_system_create(
    sir_model,
    list(
      N = N,
      I0 = I0,
      beta1 = beta1,
      beta2 = beta2,
      change_time = change_time,
      gamma = gamma
    ),
    n_particles = 1,
    deterministic = TRUE
  )

  # Run the deterministic simulation
  result <- dust2::dust_system_simulate(sys, times)

  # Extract state variables (single trajectory)
  S_vals <- result[1, 1, ]
  I_vals <- result[2, 1, ]
  R_vals <- result[3, 1, ]

  # Calculate beta at each time point for plotting
  beta_vals <- ifelse(times < change_time, beta1, beta2)

  # Create output object
  output <- list(
    times = times,
    S = S_vals,
    I = I_vals,
    R = R_vals,
    beta = beta_vals,
    parameters = list(
      N = N,
      I0 = I0,
      gamma = gamma,
      beta1 = beta1,
      beta2 = beta2,
      change_time = change_time
    )
  )

  class(output) <- c("epievolve_sir", "list")
  return(output)
}


#' Plot SIR Simulation Results
#'
#' @param x An object of class "epievolve_sir" from run_sir_step()
#' @param ... Additional arguments (not used)
#'
#' @export
plot.epievolve_sir <- function(x, ...) {

  # Set up plotting area
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))

  # Plot 1: SIR dynamics
  plot(x$times, x$S, type = "l", col = "blue", lwd = 2,
       xlab = "Time", ylab = "Number of individuals",
       main = "SIR Model Dynamics",
       ylim = c(0, max(c(x$S, x$I, x$R))))
  lines(x$times, x$I, col = "red", lwd = 2)
  lines(x$times, x$R, col = "green", lwd = 2)
  legend("right",
         legend = c("Susceptible", "Infected", "Recovered"),
         col = c("blue", "red", "green"),
         lwd = 2, bty = "n")

  # Add vertical line at change point
  abline(v = x$parameters$change_time, lty = 2, col = "gray50")

  # Plot 2: Transmission rate over time
  plot(x$times, x$beta, type = "l", col = "purple", lwd = 2,
       xlab = "Time", ylab = expression(beta),
       main = "Transmission Rate Over Time")
  abline(v = x$parameters$change_time, lty = 2, col = "gray50")

  par(mfrow = c(1, 1))
}


#' Print Summary of SIR Simulation
#'
#' @param x An object of class "epievolve_sir" from run_sir_step()
#' @param ... Additional arguments (not used)
#'
#' @export
print.epievolve_sir <- function(x, ...) {
  cat("Deterministic SIR Model with Time-Varying Transmission\n")
  cat("======================================================\n\n")

  cat("Parameters:\n")
  cat(sprintf("  Population (N):        %d\n", x$parameters$N))
  cat(sprintf("  Initial infected (I0): %d\n", x$parameters$I0))
  cat(sprintf("  Recovery rate (gamma): %.3f\n", x$parameters$gamma))
  cat(sprintf("  Beta before change:    %.3f\n", x$parameters$beta1))
  cat(sprintf("  Beta after change:     %.3f\n", x$parameters$beta2))
  cat(sprintf("  Change time:           %.1f\n\n", x$parameters$change_time))

  cat("Simulation:\n")
  cat(sprintf("  Time range:            %.1f to %.1f\n",
              min(x$times), max(x$times)))
  cat(sprintf("  Number of time points: %d\n\n", length(x$times)))

  # Final epidemic size
  final_R <- tail(x$R, 1)
  attack_rate <- final_R / x$parameters$N * 100

  cat("Epidemic outcomes:\n")
  cat(sprintf("  Final recovered:       %d\n", round(final_R)))
  cat(sprintf("  Attack rate:           %.1f%%\n", attack_rate))

  # Peak timing and size
  peak_idx <- which.max(x$I)
  peak_time <- x$times[peak_idx]
  peak_I <- x$I[peak_idx]

  cat(sprintf("  Peak infections:       %d at time %.1f\n",
              round(peak_I), peak_time))

  invisible(x)
}


#' Summary Statistics for SIR Simulation
#'
#' @param object An object of class "epievolve_sir" from run_sir_step()
#' @param ... Additional arguments (not used)
#'
#' @export
summary.epievolve_sir <- function(object, ...) {

  # Calculate R0 in each period
  R0_before <- object$parameters$beta1 / object$parameters$gamma
  R0_after <- object$parameters$beta2 / object$parameters$gamma

  # Peak statistics
  peak_idx <- which.max(object$I)
  peak_time <- object$times[peak_idx]
  peak_I <- object$I[peak_idx]

  # Final size
  final_R <- tail(object$R, 1)
  attack_rate <- final_R / object$parameters$N

  # Time to peak from start and from intervention
  time_to_peak_from_start <- peak_time
  time_to_peak_from_intervention <- peak_time - object$parameters$change_time

  result <- list(
    parameters = object$parameters,
    R0_before = R0_before,
    R0_after = R0_after,
    peak_time = peak_time,
    peak_infections = peak_I,
    time_to_peak_from_start = time_to_peak_from_start,
    time_to_peak_from_intervention = time_to_peak_from_intervention,
    final_recovered = final_R,
    attack_rate = attack_rate,
    reduction_in_transmission = (object$parameters$beta1 - object$parameters$beta2) /
      object$parameters$beta1
  )

  class(result) <- "summary.epievolve_sir"
  return(result)
}


#' Print Summary Statistics
#'
#' @param x An object of class "summary.epievolve_sir"
#' @param ... Additional arguments (not used)
#'
#' @export
print.summary.epievolve_sir <- function(x, ...) {
  cat("SIR Model Summary\n")
  cat("=================\n\n")

  cat("Basic Reproduction Numbers:\n")
  cat(sprintf("  R0 before intervention: %.2f\n", x$R0_before))
  cat(sprintf("  R0 after intervention:  %.2f\n", x$R0_after))
  cat(sprintf("  Transmission reduction: %.1f%%\n\n",
              x$reduction_in_transmission * 100))

  cat("Epidemic Dynamics:\n")
  cat(sprintf("  Peak time:              %.1f\n", x$peak_time))
  cat(sprintf("  Peak infections:        %d (%.1f%% of population)\n",
              round(x$peak_infections),
              x$peak_infections / x$parameters$N * 100))
  cat(sprintf("  Time from intervention: %.1f\n\n",
              x$time_to_peak_from_intervention))

  cat("Final Outcomes:\n")
  cat(sprintf("  Total recovered:        %d\n", round(x$final_recovered)))
  cat(sprintf("  Attack rate:            %.1f%%\n", x$attack_rate * 100))

  invisible(x)
}


# Example usage and testing
if (FALSE) {
  library(odin2)
  library(dust2)

  # Example 1: Default simulation with intervention at day 50
  result1 <- run_sir_step()
  print(result1)
  plot(result1)
  summary(result1)

  # Example 2: Strong intervention
  result2 <- run_sir_step(
    N = 10000,
    I0 = 50,
    gamma = 1/5,  # 5-day infectious period
    beta1 = 0.5,  # R0 = 2.5
    beta2 = 0.1,  # R0 = 0.5
    change_time = 30,
    times = seq(0, 200, by = 0.5)
  )
  plot(result2)

  # Example 3: Late intervention
  result3 <- run_sir_step(
    beta1 = 0.4,
    beta2 = 0.3,
    change_time = 80
  )
  plot(result3)

  # Compare different intervention timings
  par(mfrow = c(2, 2))
  for (change_t in c(20, 40, 60, 80)) {
    res <- run_sir_step(change_time = change_t, times = seq(0, 150, by = 1))
    plot(res$times, res$I, type = "l", col = "red", lwd = 2,
         main = paste("Intervention at day", change_t),
         xlab = "Time", ylab = "Infected")
    abline(v = change_t, lty = 2)
  }
}
