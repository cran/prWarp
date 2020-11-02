#' Mardia-Dryden distribution
#'
#' @description Create a matrix of 2D shape coordinates drawn from a Mardia-Dryden distribution (3D not implemented)
#'
#' @param M_ref a k x 2 refence matrix (usually the sample mean shape), where k is the number of 2D landmarks
#' @param n the number of observations
#' @param sd the standard deviation of the distribution (default = 0.02)
#'
#' @return the n x 2k matrix of shape coordinates drawn from a Mardia-Dryden distribution
#'
#' @examples
#' # 2D landmark coordinates
#' library("geomorph")
#' data("HomoMidSag")  # dataset
#' n_spec <- dim(HomoMidSag)[1]  # number of specimens
#' k <- dim(HomoMidSag)[2] / 2  # number of landmarks
#' homo_ar <- arrayspecs(HomoMidSag, k, 2)  # create an array
#'
#' # Procrustes registration
#' homo_gpa <- Morpho::procSym(homo_ar)
#' m_mshape <- homo_gpa$mshape  # average shape
#'
#' # Mardia-Dryden distribution
#' Xmd <- md.distri(m_mshape, n = n_spec, sd = 0.005)
#' # Visualization
#' plot(Xmd[, 1:k], Xmd[, (k+1):(2*k)], asp = 1, las = 1, cex = 0.5,
#' main = "Mardia-Dryden distribution", xlab = "X", ylab = "Y")
#'
#' @importFrom stats rnorm
#'
#' @export
md.distri <- function (M_ref, n, sd = 0.05) {

  if (is.null(M_ref))
    stop("supply the reference matrix 'M_ref'")
  if (!is.matrix(M_ref))
    stop("'M_ref' must be a matrix")

  if (is.null(n))
    stop("supply the number of observations 'n'")
  if (round(n) != n | n < 1)
    stop("the number of observations 'n' must be a strictly positive integer")

  if (is.null(sd))
    stop("supply the standard deviation 'sd'")
  if (sd < 0)
    stop("the standard deviation 'sd' must be a positive value")

  k <- nrow(M_ref)  # landmark number

  x_md <- matrix(NA, nrow = n, ncol = 2 * k)
  colnames(x_md) <- c(paste(1:k,"X", sep = ""), paste(1:k,"Y", sep = ""))
  for (i in 1:k) {
    x_md[, i] <- rnorm(n, mean = M_ref[i, 1], sd)  # X coordinates
    x_md[, k + i] <- rnorm(n, mean = M_ref[i, 2], sd)  # Y coordinates
  }

  return(x_md)

}
