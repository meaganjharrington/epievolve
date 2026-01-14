#' Simulate SIR Model with Two Time-Varying Transmission Rates
#'
#' @param N Total population size
#' @param I0 Initial number of infected individuals
#' @param gamma Recovery rate (1/infectious period)
#' @param beta1 Transmission rate in first period
#' @param beta2 Transmission rate in second period
#' @param beta3 Transmission rate in third period
#' @param change_time1 First change point
#' @param change_time2 Second change point
#' @param end_time End time for simulation
#'
#' @return A list containing simulation results
#'
#' @export
sir_twostep <- function(N, I0, gamma, beta1, beta2, beta3,
                        change_time1, change_time2, end_time){

  # Create time vector
  times <- seq(0, end_time, by = 1)

  # Validate inputs
  if (N <= 0) stop("N must be positive")
  if (I0 <= 0 || I0 >= N) stop("I0 must be positive and less than N")
  if (gamma <= 0) stop("gamma must be positive")
  if (beta1 <= 0 || beta2 <= 0 || beta3 <= 0) {
    stop("beta1, beta2, and beta3 must be positive")
  }
  if (change_time1 >= change_time2) {
    stop("change_time1 must be less than change_time2")
  }
  if (change_time2 > end_time) {
    stop("change_time2 must be less than or equal to end_time")
  }

  # Define the SIR model
  sir <- odin2::odin({
    deriv(S) <- -beta * S * I / N
    deriv(I) <- beta * S * I / N - gamma * I
    deriv(R) <- gamma * I

    initial(S) <- N - I0
    initial(I) <- I0
    initial(R) <- 0

    N <- parameter()
    I0 <- parameter()
    # Two step changes: if-else nested
    beta <- if (time < change_time1) beta1 else if (time < change_time2) beta2 else beta3
    beta1 <- parameter()
    beta2 <- parameter()
    beta3 <- parameter()
    change_time1 <- parameter()
    change_time2 <- parameter()
    gamma <- parameter()
  })

  pars <- list(
    N = N,
    I0 = I0,
    beta1 = beta1,
    beta2 = beta2,
    beta3 = beta3,
    change_time1 = change_time1,
    change_time2 = change_time2,
    gamma = gamma
  )

  # Create dust2 system
  sys <- dust2::dust_system_create(
    sir,
    pars,
    n_particles = 1,
    deterministic = TRUE
  )

  dust2::dust_system_set_state_initial(sys)

  # Run the deterministic simulation
  result <- dust2::dust_system_simulate(sys, times)

  # Extract state variables from matrix [state, time]
  S_vals <- result[1, ]
  I_vals <- result[2, ]
  R_vals <- result[3, ]

  # Calculate beta at each time point for plotting
  beta_vals <- ifelse(times < change_time1, beta1,
                      ifelse(times < change_time2, beta2, beta3))

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
      beta3 = beta3,
      change_time1 = change_time1,
      change_time2 = change_time2,
      end_time = end_time
    )
  )

  class(output) <- c("epievolve_sir", "list")
  return(output)
}
