#' Calculate Log-Likelihood for Marked Individuals
#'
#' This function calculates the log-likelihood for marked individuals using either a Bernoulli or Poisson observation model.
#'
#' @param obstype Character string indicating the observation type, either "bernoulli" or "poisson".
#' @param y_mark Numeric vector of observed counts for marked individuals.
#' @param K Numeric vector of sampling effort.
#' @param lambda Numeric vector of detection probabilities.
#' @param z Numeric vector of activity states (e.g., alive = 1, dead = 0).
#' @return Numeric vector of log-likelihoods.
#' @examples
#' # Example usage
#' y_mark <- c(1, 0, 2)
#' K <- c(3, 3, 3)
#' lambda <- c(0.2, 0.5, 0.1)
#' z <- c(1, 1, 0)
#' calculate_ll_bern_pois("bernoulli", y_mark, K, lambda, z)
#'
#' @export
calculate_ll_bern_pois <- function(obstype, y_mark, K, lambda, z) {
  if (obstype == "bernoulli") {
    pd <- 1 - exp(-lambda)
    ll <- dbinom(y_mark, K, pd * z, log = TRUE)
  } else {
    ll <- dpois(y_mark, K * lambda * z, log = TRUE)
  }
  return(ll)
}
