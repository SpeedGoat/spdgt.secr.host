#' Pad a Matrix with Zeros to a Specified Dimension
#'
#' This function creates a new matrix of specified dimensions,
#' filling it with zeros and placing the input matrix at a specific
#' location within the larger matrix.
#'
#' @param mat_in Input matrix to be padded
#' @param dims A numeric vector specifying the dimensions of the output matrix
#' @param position The position where the input matrix should be placed
#'        (default is the last row of the first dimension)
#'
#' @return A padded matrix with specified dimensions
#'
#' @examples
#' # Create a 3x4x5 matrix, padding a 2x4x5 input matrix at the end
#' input_mat <- array(1:40, dim = c(2, 4, 5))
#' padded_mat <- pad_matrix(input_mat, dims = c(3, 4, 5))
#'
#' # Create a matrix with input at a specific position
#' padded_mat_custom <- pad_matrix(input_mat, dims = c(3, 4, 5), position = 2)
#'
#' @export
pad_matrix <- function(mat_in, dims) {
  # Input validation
  if (missing(mat_in) || missing(dims)) {
    stop("Both 'mat_in' and 'dims' must be provided")
  }

  # Check if input matrix dimensions match the specified dims (except first dimension)
  input_dims <- dim(mat_in)
  if (length(input_dims) != length(dims)) {
    stop("Dimensions of input matrix and target dimensions must match")
  }

  # Check if other dimensions match
  if (any(input_dims[-1] != dims[-1])) {
    stop("Dimensions (except the first) must match exactly")
  }

  # Check if input matrix can fit in the target dimensions
  if (input_dims[1] > dims[1]) {
    stop("First dimension of input matrix cannot be larger than target dimensions")
  }

  # Create zero-filled matrix with target dimensions
  full_mat <- array(0, dim = dims)

  # Place input matrix in the specified position
  index_range <- 1:input_dims[1]
  full_mat[index_range, , ] <- mat_in

  return(full_mat)
}
