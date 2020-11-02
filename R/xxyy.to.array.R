#' Convert landmark data matrix into array (p x k x n)
#'
#' @description Convert a matrix of landmark coordinates into a three-dimensional array
#'
#' @param M A matrix of dimension n x [p x k] containing landmark coordinates for a set of specimens.
#' Each row contains all landmark coordinates for a single specimen.
#' The first columns correspond to the X coordinates for all landmarks, etc.
#' @param p Number of landmarks
#' @param k Number of dimensions (2 or 3)
#'
#' @return Function returns a 3D array (p x k x n),
#' where p is the number of landmark points,
#' k is the number of landmark dimensions (2 or 3),
#' and n is the number of specimens.
#' The third dimension of this array contains names for each specimen
#' if specified in the original input matrix.
#'
#' @seealso \code{\link{array.to.xxyy}}
#'
#' @examples
#' X <- matrix(rnorm(40), nrow = 4) # Random 2D coordinates of 5 landmarks for 4 specimens
#' xxyy.to.array(X, 5, 2)
#'
#' @export
xxyy.to.array <- function (M, p, k = 2) {

  if (!is.matrix(M) && !is.data.frame(M))
    stop("M must be a data frame or matrix")

  if (k < 2)
    stop("One-dimensional data cannot be used")

  if (ncol(M) != p * k)
    stop("Matrix dimensions do not match input")

  n <- nrow(M)

  # Array of dimensions p*k*n
  A <- array(NA, dim = c(p, k, n))
  for (i in 1:k) {
    A[, i, ] <- t(M[, p*(i-1) + 1:p])
  }

  dimnames(A)[[3]] <- dimnames(M)[[1]]

  return(A)

}
