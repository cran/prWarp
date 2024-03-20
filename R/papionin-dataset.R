#' papionin dataset
#'
#' 2D Cartesian coordinates of 70 landmarks
#' quantifying the skull morphology along the midsagittal plane
#' for 67 adult modern primates (mostly papionins).
#' The data correspond to a list with the 6 following elements:
#'
#' \itemize{
#' \item {\strong{coords} The 3D array of landmark coordinates}
#' \item {\strong{species} The vector of species names}
#' \item {\strong{semi_lm} The vector of semilandmark numbers of the full dataset}
#' \item {\strong{curves} The list of curves for sliding semilandmarks of the full dataset}
#' \item {\strong{links} The matrix of links between landmarks for the full dataset}
#' \item {\strong{outline} A list of 4 elements for the analysis of the outline shape: \strong{subset}, the landmark numbers for the susbet; \strong{semi_lm}, the vector of semilandmark numbers; \strong{curves}, the list of curves for sliding semilandmarks; \strong{links}, the matrix of links between landmarks for the full dataset.}
#' }
#'
#' @references Grunstra, Nicole D. S. et al. (2021)
#' Detecting phylogenetic signal and adaptation in papionin cranial shape by decomposing variation at different spatial scales.
#' \emph{Systematic Biology}, 70(4): 694--706. doi: 10.1093/sysbio/syaa093
#' @references Grunstra, Nicole D. S. et al. (2020)
#' Data form: Detecting phylogenetic signal and adaptation in papionin cranial shape by decomposing variation at different spatial scales.
#' \emph{Dryad Digital Repository}. doi: 10.5061/dryad.zkh189373
#'
#' @name papionin
#' @usage data(papionin)
#' @format A list of 6 elements.
#' @docType data
NULL
