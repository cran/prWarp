#' Principal warps
#'
#' @description Computes the principal warps and the bending energy of a reference shape configuration,
#' as well as the variance of the partial waprs, the partial warp scores and the non-affine component of shape variation
#' for 2D landmark coordinates (3D not implemented).
#' Small-scale and large-scale components of shape variation can also be computed.
#'
#' @param A a k x 2 x n array, where k is the number of 2D landmarks, and n is the sample size.
#' @param M_ref a k x 2 refence matrix (usually the sample mean shape), where k is the number of 2D landmarks
#' @param d (optional) an integer value comprised between 1 and (k-3) to compute small-scale shape components (between 1 and d)
#' and large-scale shape components (between d+1 and k-3)
#'
#' @return
#' A list containing the following named components:
#' \item{bendingEnergy}{the bending energy (the (k-3) eigenvalues of the bending energy matrix)}
#' \item{principalWarps}{the k x (k-3) matrix of principal warps (the k eigenvectors of the bending energy matrix)}
#' \item{partialWarpScores}{the n x (2k-6) matrix of partial warp (the projection of the vectors of shape coordinates, expressed as deviations from the reference shape, onto the principal warps)}
#' \item{variancePW}{the variance of the (k-3) partial warps}
#' \item{Xnonaf}{the n x 2k matrix of the non-affine component of shape variation}
#' \item{Xsmall}{the n x 2k matrix of the small-scale shape variation (if d is provided)}
#' \item{Xlarge}{the n x 2k matrix of the large-scale shape variation (if d is provided)}
#'
#' @seealso See \code{\link[Morpho]{CreateL}} for the creation of the bending energy matrix
#'
#' @references Bookstein FL. (1989).
#' Principal Warps: Thin-plate splines and the decomposition of deformations.
#' \emph{IEEE Transactions on pattern analysis and machine intelligence 11(6)}: 567--585.
#' \url{https://ieeexplore.ieee.org/abstract/document/24792}
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
#' m_overall <- homo_gpa$rotated  # Procrustes coordinates
#' m_mshape <- homo_gpa$mshape  # average shape
#'
#' # Computation of bending energy, partial warp scores etc.
#' homo_be_pw <- create.pw.be(m_overall, m_mshape)
#'
#' # Partial warp variance as a function of bending energy
#' logInvBE <- log((homo_be_pw$bendingEnergy)^(-1))  # inverse log bending energy
#' logPWvar <- log(homo_be_pw$variancePW)  # log variance of partial warps
#' mod <- lm(logPWvar ~ logInvBE)  # linear regression
#' # Plot log PW variance on log BE^-1 with regression line
#' plot(logInvBE, logPWvar, col = "white", asp = 1,
#' main = "PW variance against inverse BE", xlab = "log 1/BE", ylab = "log PW variance")
#' text(logInvBE, logPWvar, labels = names(logPWvar), cex = 0.5)
#' abline(mod, col = "blue")
#'
#' @importFrom Morpho CreateL
#' @importFrom stats var
#'
#' @export
create.pw.be <- function (A, M_ref, d = NULL) {

  if (is.null(A))
    stop("supply the array of landmark coordinates 'A'")
  if (!is.array(A))
    stop("'A' must be an array")
  if (length(dim(A)) != 3)
    stop("'A' must be a three-dimensional array")
  if (dim(A)[2] != 2)
    stop("'A' must be a k x 2 x n array")

  if (is.null(M_ref))
    stop("supply the reference matrix 'M_ref'")
  if (!is.matrix(M_ref))
    stop("'M_ref' must be a matrix")
  if (any(dim(A)[1:2] != dim(M_ref)))
    stop("A and M_ref must have the same number of 2D landmarks")

  k <- dim(A)[[1]]  # landmark number
  nspec <- dim(A)[[3]]  # specimen number

  # Create the bending energy matrix as specified in Bookstein (1989)
  m_L <- CreateL(M_ref)
  m_kxk <- as.matrix(m_L$Lsubk)  # bending energy matrix

  # Take only the (k-3) non-zero eigenvalues and the corresponding eigenvectors
  m_eigen <- eigen(m_kxk)
  m_PW <- m_eigen$vectors[, 1:(k-3)]  # principal warps
  m_be <- m_eigen$values[1:(k-3)]  # bending energy

  # Deviation of the landmark coordinates from the reference shape
  m_x <- matrix(NA, nrow = nspec, ncol = k)
  m_y <- matrix(NA, nrow = nspec, ncol = k)
  for (i in 1:nspec){
    m_x[i, ] <- A[, 1, i] - M_ref[, 1]
    m_y[i, ] <- A[, 2, i] - M_ref[, 2]
  }

  # Partial warps
  m_x_PWscores <- m_x %*% m_PW  # partial warp scores for the X dimension
  m_y_PWscores <- m_y %*% m_PW  # partial warp scores for the Y dimension
  colnames(m_x_PWscores) <- paste(1:(k-3), "X", sep = "")
  colnames(m_y_PWscores) <- paste(1:(k-3), "Y", sep = "")
  m_PWscores <- cbind(m_x_PWscores, m_y_PWscores)
  if (!is.null(dimnames(A)[[3]])) {
    rownames(m_PWscores) <- dimnames(A)[[3]]
  }

  # Non-affine component of shape variation
  m_x_nonaf <- m_x_PWscores %*% t(m_PW)
  m_y_nonaf <- m_y_PWscores %*% t(m_PW)
  for (i in 1:nspec){
    m_x_nonaf[i, ] <- m_x_nonaf[i, ] + M_ref[, 1] # add the reference shape
    m_y_nonaf[i, ] <- m_y_nonaf[i, ] + M_ref[, 2] # add the reference shape
  }
  m_nonaf <- cbind(m_x_nonaf, m_y_nonaf)

  # Variance of each partial warp (mean of the variances for X and Y coordinates)
  varPW <- rep(NA, (k-3))
  for (j in 1:(k-3)) {
    var_x <- var(m_PWscores[, j])
    var_y <- var(m_PWscores[, (k - 3 + j)])
    varPW[j] <- (var_x + var_y) / 2
  }
  names(varPW) <- paste("PW", 1:(k-3), sep = "")

  # Small- and large-scale shape components
  m_small <- NULL
  m_large <- NULL
  if (!is.null(d)) {
    if (round(d) != d) {
      warning("d must be an integer for the computation of small- and large-scale shape components.")
    }
    if (round(d) == d & (d < 2 | d > k - 4)) {
      warning("d must be comprised between 2 and (k-4) for the computation of small- and large-scale shape components.")
    }
    if (round(d) == d & d > 1 & d < k - 3) {
      # Small-scale component of shape variation
      m_x_s <- m_x_PWscores[,1:d] %*% t(m_PW[,1:d])
      m_y_s <- m_y_PWscores[,1:d] %*% t(m_PW[,1:d])
      for (i in 1:nspec){
        m_x_s[i, ] <- m_x_s[i, ] + M_ref[, 1] # add the reference shape
        m_y_s[i, ] <- m_y_s[i, ] + M_ref[, 2] # add the reference shape
      }
      m_small <- cbind(m_x_s, m_y_s)
      # Large-scale component of shape variation
      m_x_l <- m_x_PWscores[,(d+1):(k-3)] %*% t(m_PW[,(d+1):(k-3)])
      m_y_l <- m_y_PWscores[,(d+1):(k-3)] %*% t(m_PW[,(d+1):(k-3)])
      for (i in 1:nspec){
        m_x_l[i, ] <- m_x_l[i, ] + M_ref[, 1] # add the reference shape
        m_y_l[i, ] <- m_y_l[i, ] + M_ref[, 2] # add the reference shape
      }
      m_large <- cbind(m_x_l, m_y_l)
    }
  }

  # Output
  results <- list(bendingEnergy = m_be,
                  principalWarps = m_PW,
                  partialWarpScores = m_PWscores,
                  variancePW = varPW,
                  Xnonaf = m_nonaf,
                  Xsmall = m_small,
                  Xlarge = m_large)

  return(results)

}
