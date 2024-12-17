#' Check if a Point is Inside a Study Area
#'
#' Determines whether a given point is within the boundaries of a study area,
#' either using rectangular limits or a more complex polygon.
#'
#' @param point A numeric vector of length 2 representing (x, y) coordinates
#' @param xlim A numeric vector of length 2 specifying the x-axis limits (min, max)
#' @param ylim A numeric vector of length 2 specifying the y-axis limits (min, max)
#' @param vertices Optional polygon vertices for complex boundary checking
#' @param use_vertices Logical. If TRUE, uses polygon vertices for boundary check.
#'        If FALSE, uses rectangular xlim and ylim.
#'
#' @return Logical. TRUE if the point is inside the specified area, FALSE otherwise.
#'
#' @examples
#' # Rectangular boundary check
#' point_in_area(c(5, 5), xlim = c(0, 10), ylim = c(0, 10))
#'
#' # Polygon boundary check (requires sp package)
#' library(sp)
#' poly_vertices <- cbind(x = c(0, 10, 10, 0), y = c(0, 0, 10, 10))
#' point_in_area(c(5, 5), vertices = poly_vertices, use_vertices = TRUE)
#'
#' @export
point_in_area <- function(point,
                          xlim = NULL,
                          ylim = NULL,
                          vertices = NULL,
                          use_vertices = FALSE) {
  # Validate inputs
  if (!is.numeric(point) || length(point) != 2) {
    stop("point must be a numeric vector of length 2")
  }

  # Check using rectangular boundaries
  if (!use_vertices) {
    # Require both xlim and ylim for rectangular check
    if (is.null(xlim) || is.null(ylim) ||
        length(xlim) != 2 || length(ylim) != 2) {
      stop("For rectangular check, provide valid xlim and ylim")
    }

    return(
      point[1] < xlim[2] &
        point[1] > xlim[1] &
        point[2] < ylim[2] &
        point[2] > ylim[1]
    )
  }

  # Check using polygon vertices
  if (use_vertices) {
    # Require vertices
    if (is.null(vertices)) {
      stop("When use_vertices = TRUE, vertices must be provided")
    }

    # Recommend sp package if not already in namespace
    if (!requireNamespace("sp", quietly = TRUE)) {
      warning("sp package recommended for polygon point-in-polygon checks")
    }

    # Determine if point is in polygon
    return(in_poly(point, vertices))
  }
}

#' Determine if a Point is Inside a Polygon
#'
#' Checks whether a given point is inside a polygon using the ray-casting algorithm.
#'
#' @param point A numeric vector of length 2 representing (x, y) coordinates of the point to check
#' @param vertices A numeric matrix or data frame with two columns (x and y) representing polygon vertices
#'
#' @return Logical. TRUE if the point is inside the polygon, FALSE otherwise.
#'
#' @details
#' This function implements the ray-casting algorithm to determine if a point is inside a polygon.
#' The algorithm works by casting a horizontal ray from the point and counting the number of
#' times it intersects the polygon edges. An odd number of intersections indicates the point
#' is inside the polygon.
#'
#' @examples
#' # Define a simple triangle
#' vertices <- matrix(c(0,0, 10,0, 5,10), ncol=2, byrow=TRUE)
#'
#' # Point inside the polygon
#' point_inside <- c(5, 5)
#' inout(point_inside, vertices)  # Should return TRUE
#'
#' # Point outside the polygon
#' point_outside <- c(20, 20)
#' inout(point_outside, vertices)  # Should return FALSE
#'
#' @importFrom stats pnorm
in_poly <- function(point, vertices) {
  # Input validation
  if (!is.numeric(point) || length(point) != 2) {
    stop("point must be a numeric vector of length 2")
  }

  if (!is.matrix(vertices) && !is.data.frame(vertices)) {
    stop("vertices must be a matrix or data frame")
  }

  if (ncol(vertices) != 2) {
    stop("vertices must have exactly two columns (x and y)")
  }

  # Ensure vertices form a closed polygon by appending the first vertex to the end
  vertices <- rbind(vertices, vertices[1,])

  # Ray-casting algorithm
  intersect_count <- 0
  x <- point[1]
  y <- point[2]

  for (i in 1:(nrow(vertices) - 1)) {
    # Get current edge vertices
    x1 <- vertices[i, 1]
    y1 <- vertices[i, 2]
    x2 <- vertices[i+1, 1]
    y2 <- vertices[i+1, 2]

    # Check if the ray intersects the edge
    # Condition 1: y-coordinate is within the edge's y range
    # Condition 2: horizontal ray from point crosses the edge
    if (((y1 > y) != (y2 > y)) &&
        (x < (x2 - x1) * (y - y1) / (y2 - y1) + x1)) {
      intersect_count <- intersect_count + 1
    }
  }

  # Odd number of intersections means point is inside
  return(intersect_count %% 2 == 1)
}
