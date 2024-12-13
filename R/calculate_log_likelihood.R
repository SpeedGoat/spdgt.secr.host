#' Calculate Log-Likelihood
#'
#' Computes the log-likelihood given a value, center, bounds, and sigma.
#' @param value The observed value.
#' @param center The center of the distribution.
#' @param bounds A numeric vector of length 2 specifying bounds.
#' @param sigma The standard deviation.
#' @return The log-likelihood.
#'
#' @examples
#' # Example inputs
#' value <- 1.5               # Observed value
#' center <- 1.0              # Mean of the distribution
#' bounds <- c(0, 2)          # Bounds of the truncated normal distribution
#' sigma <- 0.5               # Standard deviation
#'
#' # Call the function
#' log_likelihood <- calculate_log_likelihood(value, center, bounds, sigma)
#'
#' # Print the result
#' print(paste("Log-likelihood for value =", value, "is:", log_likelihood))
#'
#' @export
calculate_log_likelihood <- function(value, center, bounds, sigma) {
  log(
    dnorm(value, center, sigma) /
      (pnorm(bounds[2], center, sigma) - pnorm(bounds[1], center, sigma))
  )
}

