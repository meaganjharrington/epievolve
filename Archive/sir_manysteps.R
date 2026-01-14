#' Simulate SIR Model with Multiple Time-Varying Transmission Rates
#'
#' @param N Total population size
#' @param I0 Initial number of infected individuals
#' @param gamma Recovery rate (1/infectious period)
#' @param time_points Vector of time points where beta changes
#' @param beta_values Vector of beta values corresponding to time_points
#' @param end_time End time for simulation
#' @param time_step Time step for output (default: 0.25)
#'
#' @return A list containing:
#'   \item{times}{Time points}
#'   \item{S}{Susceptible counts over time}
#'   \item{I}{Infected counts over time}
#'   \item{R}{Recovered counts over time}
#'   \item{beta}{Transmission rate over time}
#'   \item{parameters}{List of parameters used}
#'
#' @export
sir_manysteps <- function(N, I0, gamma, time_points, beta_values,
                          end_time, time_step = 0.25){

  # Create time vector
  times <- seq(0, end_time, by = time_step)

  # Validate inputs
  if (N <= 0) stop("N must be positive")
  if (I0 <= 0 || I0 >= N) stop("I0 must be positive and less than N")
  if (gamma <= 0) stop("gamma must be positive")
  if (length(time_points) != length(beta_values)) {
    stop("time_points and beta_values must have the same length")
  }
  if (any(beta_values < 0)) stop("All beta_values must be non-negative")
  if (!all(diff(time_points) > 0)) {
    stop("time_points must be strictly increasing")
  }
  if (min(time_points) > 0) {
    stop("time_points must start at 0 or earlier")
  }
  if (max(time_points) < end_time) {
    warning("Last time_point is before end_time; beta will be constant after last time_point")
  }

  # Define the SIR model
  sir <- odin2::odin({
    deriv(S) <- -beta * S * I / N
    deriv(I) <- beta * S * I / N - gamma * I
    deriv(R) <- gamma * I

    initial(S) <- N - I0
    initial(I) <- I0
    initial(R) <- 0

    N     <- parameter()
    I0    <- parameter()
    gamma <- parameter()

    # Define arrays for interpolation
    dim(time_points) <- parameter(rank = 1)
    dim(beta_values) <- length(time_points)

    beta <- interpolate(time_points, beta_values, "constant")
  })

  # Create parameter list
  pars <- list(
    N = N,
    I0 = I0,
    gamma = gamma,
    time_points = time_points,
    beta_values = beta_values
  )

  # Create dust2 system
  sys <- dust2::dust_system_create(
    sir,
    pars,
    n_particles = 1,
    deterministic = TRUE
  )

  # Set initial state
  dust2::dust_system_set_state_initial(sys)

  # Run the deterministic simulation
  result <- dust2::dust_system_simulate(sys, times)

  # Extract state variables from matrix [state, time]
  S_vals <- result[1, ]
  I_vals <- result[2, ]
  R_vals <- result[3, ]

  # Interpolate beta values for plotting
  # Use constant interpolation (step function)
  beta_vals <- approx(
    x = time_points,
    y = beta_values,
    xout = times,
    method = "constant",
    rule = 2,  # Use nearest value outside range
    f = 0      # Right-continuous step function
  )$y

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
      time_points = time_points,
      beta_values = beta_values,
      end_time = end_time,
      time_step = time_step
    )
  )

  class(output) <- c("epievolve_sir", "list")
  return(output)
}
