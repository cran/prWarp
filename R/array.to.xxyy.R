#' Convert (p x k x n) data array into 2D data matrix
#'
#' @description Convert a three-dimensional array of landmark coordinates into a two-dimensional matrix
#'
#' @param A A 3D array (p x k x n) containing landmark coordinates for a set of specimens
#'
#' @return Function returns a two-dimensional matrix of dimension (n x [p x k]),
#' where rows represent specimens and columns represent variables.
#' The p first columns correspond to X coordinates, etc.
#'
#' @seealso \code{\link{xxyy.to.array}}
#'
#' @examples
#' A <- array(rnorm(40), c(5, 2, 4))  # random 2D coordinates of 5 landmarks for 4 specimens
#' array.to.xxyy(A)
#'
#' @export
array.to.xxyy <- function (A) {

  if (!is.array(A) | length(dim(A)) != 3)
    stop("'A' must be a 3D array (see 'array.to.xxyy').")

  M <- t(A[, 1, ])  # X coordinates
  for (i in 2:dim(A)[2]) {
    M <- cbind(M, t(A[, i, ]))
  }

  rownames(M) <- dimnames(A)[[3]]

  cnames <- c(paste(1:dim(A)[1], "X", sep = ""),
              paste(1:dim(A)[1], "Y", sep = ""))

  if (dim(A)[2] == 2)
    colnames(M) <- cnames
  if (dim(A)[2] == 3)
    colnames(M) <- c(cnames, paste(1:dim(A)[1], "Z", sep = ""))

  return(M)

}
