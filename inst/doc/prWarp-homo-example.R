## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, out.width = '100%', dpi = 600)

## ----packages-----------------------------------------------------------------
library("prWarp")

## ----data---------------------------------------------------------------------
data("HomoMidSag")
k <- dim(HomoMidSag)[2] / 2  # number of landmarks
n_spec <- dim(HomoMidSag)[1]  # number of specimens
homo_ar <- geomorph::arrayspecs(HomoMidSag, k, 2)  # create an array
dimnames(homo_ar)[[1]] <- 1:k
dimnames(homo_ar)[[2]] <- c("X", "Y")

## ----gpa, include=FALSE-------------------------------------------------------
homo_gpa <- Morpho::procSym(homo_ar)
m_overall <- homo_gpa$rotated  # Procrustes coordinates
m_mshape <- homo_gpa$mshape  # average shape

## ----plot mhape---------------------------------------------------------------
plot(m_mshape, asp = 1, main = "Average shape", xlab = "X", ylab = "Y")

## ----partial warp decomposition-----------------------------------------------
homo_be_pw <- create.pw.be(m_overall, m_mshape)

## ----PW variance vs BE--------------------------------------------------------
# Computation of log BE^-1 for the (k-3) partial warps
logInvBE <- log((homo_be_pw$bendingEnergy)^(-1))
# Computation of log PW variance for the (k-3) partial warps
logPWvar <- log(homo_be_pw$variancePW)
# Linear regression of the log PW variance on the log BE^-1
mod <- lm(logPWvar ~ logInvBE)
# Plot log PW variance on log BE^-1 with regression line
plot(logInvBE, logPWvar, col = "white", asp = 1, main = "PW variance against inverse BE", sub = paste("slope =", round(mod$coefficients[2], 2)), xlab = "log 1/BE", ylab = "log PW variance")
text(logInvBE, logPWvar, labels = names(logPWvar), cex = 0.5)
abline(mod, col = "blue")

## ----plot non aff-------------------------------------------------------------
# Compute the trace of t(Xnonaf) %*% Xnonaf
tr_nonaf <- sum(diag(t(homo_be_pw$Xnonaf) %*% homo_be_pw$Xnonaf))
# Convert matrix into a 3D array
Anonaf <- xxyy.to.array(homo_be_pw$Xnonaf, p = k, k = 2) 
# Plot the non-affine shape variation around the mean
geomorph::plotAllSpecimens(Anonaf, plot.param = list(pt.cex = 0.3, mean.cex = 0.8, mean.col = "red"))

## -----------------------------------------------------------------------------
# Compute the self-similar distribution
Xdefl <- ssim.distri(m_mshape, n = n_spec, sd = 0.05, f = 1)
# Compute the trace of t(Xdefl) %*% Xdefl
tr_defl <- sum(diag(t(Xdefl) %*% Xdefl))
# Convert matrix into a 3D array
Adefl <- xxyy.to.array(Xdefl, p = k, k = 2) 
# Plot the self-similar distribution
geomorph::plotAllSpecimens(Adefl, plot.param = list(pt.cex = 0.3, mean.cex = 0.8, mean.col = "red"))


## -----------------------------------------------------------------------------
# Compute the Mardia-Dryden distribution
Xmd <- md.distri(m_mshape, n = n_spec, sd = 0.005)
# Convert matrix into a 3D array
Amd <- xxyy.to.array(Xmd, p = k, k = 2) 
# Plot the Mardia-Dryden distribution
geomorph::plotAllSpecimens(Amd, plot.param = list(pt.cex = 0.3, mean.cex = 0.8, mean.col = "red"))

