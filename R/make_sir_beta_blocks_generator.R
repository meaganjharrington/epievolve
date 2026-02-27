make_sir_beta_blocks_generator <- function() {
  odin2::odin({
    # initial
    initial(S) <- N - I0
    initial(I) <- I0
    initial(incidence, zero_every = 1) <- 0

    # deterministic transition probabilities
    p_SI <- 1 - exp(-beta_t * I / N * dt)
    p_IR <- 1 - exp(-gamma * dt)

    new_inf <- S * p_SI
    new_rec <- I * p_IR

    update(incidence) <- new_inf
    update(S) <- S - new_inf
    update(I) <- I + new_inf - new_rec

    # parameters
    N <- parameter()
    I0 <- parameter()
    gamma <- parameter()
    n_beta <- parameter()

    beta_times  <- parameter()
    beta_values <- parameter()
    dim(beta_times) <- n_beta
    dim(beta_values) <- n_beta

    beta_t <- interpolate(beta_times, beta_values, "constant")

    # data likelihood
    cases <- data()
    cases ~ Poisson(incidence)
  })
}
