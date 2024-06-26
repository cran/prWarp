% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tps.all.R
\name{tps.all}
\alias{tps.all}
\title{Thin plate spline mapping (2D and 3D) for several sets of landmark coordinates}
\usage{
tps.all(X_array, REF_array, TAR_matrix)
}
\arguments{
\item{X_array}{original coordinates
- a 3D array (p x k x n) containing original landmark coordinates for a set of specimens}

\item{REF_array}{reference coordinates (e.g., outline landmarks for all specimens)
- a 3D array (p x k x n) containing reference landmark coordinates for a set of specimens}

\item{TAR_matrix}{target coordinates (e.g., average outline landmarks)
- a matrix (n x k) containing target landmark coordinates}
}
\value{
Function returns a 3D array (p x k x n) containing the deformed input
(original landmark set warped onto the target matrix).
}
\description{
Maps landmarks via thin plate spline
based on a reference and a target configuration in 2D or in 3D.
This function is an extension of the tps3d function for a set of specimens.
}
\details{
p is the number of landmark points, k is the number of landmark dimensions (2 or 3), and n is the number of specimens.
}
\examples{
data("papionin")  # load dataset
# Full dataset: 70 landmarks
papionin_ar <- papionin$coords
# Outline dataset: subset of 54 landmarks
outline_ar <- papionin_ar[papionin$outline$subset, , ]

# Subset: Macaca only
mac <- grep("Macaca", papionin$species)  # genus Macaca
papionin_macaca <- papionin_ar[, , mac]
outline_macaca <- outline_ar[, , mac]

# Landmark sliding by minimizing bending energy + superimposition (GPA)
library("Morpho")
papionin_gpa <- procSym(papionin_macaca, SMvector = papionin$semi_lm,
                        outlines = papionin$curves)
outline_gpa <- procSym(outline_macaca, SMvector = papionin$outline$semi_lm,
                       outlines = papionin$outline$curves)

# Warping the slid landmarks of the full landmark dataset to the average outline shape
residual_shape <- tps.all(X_array = papionin_gpa$dataslide,
                          REF_array = outline_gpa$dataslide,
                          TAR_matrix = outline_gpa$mshape)

}
\references{
Bookstein FL. (1989).
Principal Warps: Thin-plate splines and the decomposition of deformations.
\emph{IEEE Transactions on pattern analysis and machine intelligence}, 11(6): 567--585.
\url{https://ieeexplore.ieee.org/abstract/document/24792}
}
\seealso{
See \code{\link[Morpho]{tps3d}}
}
