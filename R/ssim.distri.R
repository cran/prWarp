#' Self-similar distribution
#'
#' @description Create a matrix of 2D shape coordinates drawn from a self-similar distribution (3D not implemented)
#'
#' @param M_ref a k x 2 refence matrix (usually the sample mean shape), where k is the number of 2D landmarks
#' @param n the number of observations
#' @param sd the standard deviation of the distribution (default = 0.02)
#' @param f a scaling factor
#'
#' @return the n x 2k matrix of shape coordinates drawn from a self-similar distribution
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
#' # Self-similar distribution
#' Xdefl <- ssim.distri(m_mshape, n = n_spec, sd = 0.05, f = 1)
#' # Visualization
#' plot(Xdefl[, 1:k], Xdefl[, (k+1):(2*k)], asp = 1, las = 1, cex = 0.5,
#' main = "Self-similar distribution", xlab = "X", ylab = "Y")
#'
#' @importFrom stats rnorm
#'
#' @export
ssim.distri <- function (M_ref, n, sd = 0.02, f = 1) {

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

  if (is.null(f))
    stop("supply the scaling factor 'f'")
  if (f < 0)
    stop("the scaling factor 'sd' must be a positive value")

  k <- nrow(M_ref)  # landmark number

  # Create the bending energy matrix as specified in Bookstein (1989)
  m_L <- CreateL(M_ref)
  m_kxk <- as.matrix(m_L$Lsubk)  # bending energy (BE) matrix

  # Take only the (k-3) non-zero eigenvalues and the corresponding eigenvectors
  m_eigen <- eigen(m_kxk)
  m_PW <- m_eigen$vectors[, 1:(k-3)]  # eigenvectors of the BE matrix
  m_be <- m_eigen$values[1:(k-3)]  # eigenvalues of the BE matrix

  # Compute a Mardia-Dryden distribution and remove reference shape
  m_x <- matrix(NA, nrow = n, ncol = k)
  m_y <- matrix(NA, nrow = n, ncol = k)
  for (i in 1:k) {
    m_x[, i] <- rnorm(n, mean = 0, sd)
    m_y[, i] <- rnorm(n, mean = 0, sd)
  }

  # Self-similar distribution
  m_x_defl <- m_x %*% m_PW %*% diag(m_be ^ (-0.5)) %*% t(m_PW)
  m_y_defl <- m_y %*% m_PW %*% diag(m_be ^ (-0.5)) %*% t(m_PW)
  for (i in 1:n){
    m_x_defl[i, ] <- f * m_x_defl[i, ] + M_ref[, 1] # add the reference shape
    m_y_defl[i, ] <- f * m_y_defl[i, ] + M_ref[, 2] # add the reference shape
  }
  colnames(m_x_defl) <- paste(1:k, "X", sep = "")
  colnames(m_y_defl) <- paste(1:k, "Y", sep = "")
  Xdefl <- cbind(m_x_defl, m_y_defl)

  return(Xdefl)

}
