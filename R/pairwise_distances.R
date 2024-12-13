#' Calculate Pairwise Distances Between 2D Points
#'
#' This function computes pairwise Euclidean distances between rows of two 2D matrices,
#' data frames, or vectors treated as single-row matrices.
#'
#' @param x A 2-column matrix, data frame, or a length-2 vector representing the first set of points.
#' @param y Optional. A 2-column matrix, data frame, or a length-2 vector representing the second set of points.
#' If NULL, distances are calculated within `x`.
#' @return A matrix of pairwise distances where rows correspond to points in `x` and columns correspond to points in `y`.
#' @examples
#' # Example with matrices
#' x <- matrix(c(1, 2, 3, 4), ncol = 2)
#' y <- matrix(c(5, 6, 7, 8), ncol = 2)
#' pairwise_distances(x, y)
#'
#' # Example with a vector
#' x <- c(1, 2)
#' pairwise_distances(x)
#'
#' @export
pairwise_distances <- function(x, y = NULL) {
  # Ensure x is a 2-column matrix
  if (is.null(dim(x)) && length(x) == 2) {
    x <- matrix(x, nrow = 1)
  }
  if (ncol(x) != 2) {
    stop("Argument 'x' must be a 2-column matrix or data frame, or a length-2 vector.", call. = FALSE)
  }

  # If y is NULL, set y to x
  if (is.null(y)) {
    y <- x
  } else {
    # Ensure y is a 2-column matrix
    if (is.null(dim(y)) && length(y) == 2) {
      y <- matrix(y, nrow = 1)
    }
    if (ncol(y) != 2) {
      stop("Argument 'y' must be a 2-column matrix or data frame, or a length-2 vector.", call. = FALSE)
    }
  }

  # Compute pairwise distances
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = FALSE)
}
