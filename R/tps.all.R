#' Thin plate spline mapping (2D and 3D) for several sets of landmark coordinates
#'
#' @description Maps landmarks via thin plate spline
#' based on a reference and a target configuration in 2D and 3D.
#' This function is an extension of the tps3d function a set of specimens.
#'
#' @param X_array original coordinates
#' - a 3D array (p x k x n) containing original landmark coordinates for a set of specimens
#' @param REF_array reference coordinates (e.g., outline landmarks for all specimens)
#' - a 3D array (p x k x n) containing reference landmark coordinates for a set of specimens
#' @param TAR_matrix target coordinates (e.g., average outline landmarks)
#' - a matrix (n x k) containing target landmark coordinates
#'
#' @details p is the number of landmark points, k is the number of landmark dimensions (2 or 3), and n is the number of specimens.
#'
#' @seealso See \code{\link[Morpho]{tps3d}}
#'
#' @references Bookstein FL. (1989).
#' Principal Warps: Thin-plate splines and the decomposition of deformations.
#' \emph{IEEE Transactions on pattern analysis and machine intelligence 11(6)}: 567--585.
#' \url{https://ieeexplore.ieee.org/abstract/document/24792}
#'
#' @return Function returns a 3D array (p x k x n) containing the deformed input
#' (original landmark set warped onto the target matrix).
#'
#' @examples
#' data("papionin")  # load dataset
#' # Full dataset: 70 landmarks
#' papionin_ar <- papionin$coords
#' # Outline dataset: subset of 54 landmarks
#' outline_ar <- papionin_ar[papionin$outline$subset, , ]
#'
#' # Subset: Macaca only
#' mac <- grep("Macaca", papionin$species)  # genus Macaca
#' papionin_macaca <- papionin_ar[, , mac]
#' outline_macaca <- outline_ar[, , mac]
#'
#' # Landmark sliding by minimizing bending energy + superimposition (GPA)
#' library("Morpho")
#' papionin_gpa <- procSym(papionin_macaca, SMvector = papionin$semi_lm,
#'                         outlines = papionin$curves)
#' outline_gpa <- procSym(outline_macaca, SMvector = papionin$outline$semi_lm,
#'                        outlines = papionin$outline$curves)
#'
#' # Warping the slid landmarks of the full landmark dataset to the average outline shape
#' residual_shape <- tps.all(X_array = papionin_gpa$dataslide,
#'                           REF_array = outline_gpa$dataslide,
#'                           TAR_matrix = outline_gpa$mshape)
#'
#' @importFrom Morpho tps3d
#'
#' @export
tps.all <- function (X_array, REF_array, TAR_matrix) {

  if (!is.array(X_array) | length(dim(X_array)) != 3)
    stop("'X_array' must be a 3D array.")
  if (!is.array(REF_array) | length(dim(REF_array)) != 3)
    stop("'REF_array' must be a 3D array.")
  if (!is.matrix(TAR_matrix))
    stop("'TAR_matrix' must be a matrix.")
  if (dim(X_array)[3] != dim(REF_array)[3])
    stop("'X_array' and 'REF_array' must have the same number of dimensions.")
  if (any(dim(REF_array)[1:2] != dim(TAR_matrix)))
    stop("Dimensions of landmark data in 'TAR_matrix' do not correspond to the reference array.")

  X_warp <- array(NA, dim = dim(X_array), dimnames = dimnames(X_array))
  for (i in 1:dim(X_array)[3]) {
    X_warp[, , i] <- tps3d(X_array[,,i], REF_array[,,i], TAR_matrix)
  }

  return(X_warp)

}
