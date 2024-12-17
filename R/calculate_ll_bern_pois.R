#' Calculate Log-Likelihood
#'
#' This function calculates the log-likelihood using either a Bernoulli or Poisson observation model.
#'
#' @param obstype Character string indicating the observation type, either "bernoulli" or "poisson".
#' @param y Numeric vector of observed counts.
#' @param K Numeric vector of sampling effort.
#' @param lambda Numeric vector of detection probabilities.
#' @param z Numeric vector of activity states (e.g., alive = 1, dead = 0).
#' @return Numeric vector of log-likelihoods.
#' @examples
#' # Example usage
#' y <- c(1, 0, 2)
#' K <- c(3, 3, 3)
#' lambda <- c(0.2, 0.5, 0.1)
#' z <- c(1, 1, 0)
#' calculate_ll_bern_pois("bernoulli", y, K, lambda, z)
#'
#' @export
calculate_ll_bern_pois <- function(obstype, y, K, lambda, z) {
  if (obstype == "bernoulli") {
    pd <- 1 - exp(-lambda)
    ll <- dbinom(y, K, pd * z, log = TRUE)
  } else {
    ll <- dpois(y, K * lambda * z, log = TRUE)
  }
  return(ll)
}
