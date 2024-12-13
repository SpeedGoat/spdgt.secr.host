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
#' @export
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
