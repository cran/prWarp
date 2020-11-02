## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, out.width = '100%', dpi = 600)

## ----packages, message = FALSE------------------------------------------------
library("prWarp")
library("geomorph")

## ----data---------------------------------------------------------------------
data("papionin")  # load dataset
species <- papionin$species  # species
papionin_ar <- papionin$coords  # landmark coordinates
k <- dim(papionin_ar)[1]  # number of landmarks
n_spec <- dim(papionin_ar)[3]  # number of specimens
n_species <- length(levels(species))  # number of species

## ----all semilandmarks--------------------------------------------------------
# Semilandmarks
all_semi_lm <- papionin$semi_lm
# Curves for the sliding process
all_curves <- papionin$curves
# Links between landmarks for visualization
all_links <- papionin$links

## ----plot all fixed vs sliding landmarks--------------------------------------
plot(papionin_ar[,,1], pch = 16, col = "black", cex = 0.7, asp = 1)
points(papionin_ar[-all_semi_lm,,1], pch = 16, col = "red")
text(papionin_ar[-all_semi_lm,,1], label = c(1:k)[-all_semi_lm], pos = 1, col = "red", cex = 0.8)
text(papionin_ar[all_semi_lm,,1], label = all_semi_lm, pos = 1, col = "black", cex = 0.6)

## ----all gpa sliding, message = FALSE-----------------------------------------
papionin_gpa <- Morpho::procSym(papionin_ar, SMvector = all_semi_lm, outlines = all_curves, pcAlign = FALSE)  # sliding + GPA
papionin_ar_slid <- papionin_gpa$dataslide  # slid landmarks in the original space
total_shape <- papionin_gpa$rotated  # Procrustes coordinates
total_shape_average <- papionin_gpa$mshape  # average shape

## ----plot all fixed vs slid landmarks-----------------------------------------
plot(papionin_ar_slid[,,1], pch = 16, col = "black", cex = 0.7, asp = 1)
points(papionin_ar_slid[-all_semi_lm,,1], pch = 16, col = "red")
text(papionin_ar_slid[-all_semi_lm,,1], label = c(1:k)[-all_semi_lm], pos = 1, col = "red", cex = 0.8)
text(papionin_ar_slid[all_semi_lm,,1], label = all_semi_lm, pos = 1, col = "black", cex = 0.6)

## ----plot total shape---------------------------------------------------------
plotAllSpecimens(total_shape, mean = FALSE, plot.param = list(pt.cex = 0.5))
for (i in 1:n_spec) {
  segments(total_shape[all_links[,1],1,i], 
           total_shape[all_links[,1],2,i], 
           total_shape[all_links[,2],1,i], 
           total_shape[all_links[,2],2,i])
}

## ----total shape between species----------------------------------------------
# Computation of the average total shape for each species
total_shape_between <- array(NA, dim = c(k, 2, n_species))
for (i in 1:n_species) {
  spi <- which(species == levels(species)[i])
  total_shape_between[,,i] <- mshape(total_shape[, , spi])
}
# Plot total shape between species
plotAllSpecimens(total_shape_between, mean = FALSE, plot.param = list(pt.cex = 0.5))
for (i in 1:n_species) {
  segments(total_shape_between[all_links[,1],1,i], 
           total_shape_between[all_links[,1],2,i], 
           total_shape_between[all_links[,2],1,i], 
           total_shape_between[all_links[,2],2,i])
}

## ----total shape within species-----------------------------------------------
# Computation of the pooled individual within-species total shape
total_shape_within <- array(NA, dim = c(k, 2, n_spec))
for (i in 1:n_species) {
  spi <- which(species == levels(species)[i])
  for (j in 1:length(spi)) {
    total_shape_within[,,spi[j]] <- total_shape[, , spi[j]] - total_shape_between[,,i] + mshape(total_shape_between)
  }
}
# Plot the pooled individual within-species total shape
plotAllSpecimens(total_shape_within, mean = FALSE, plot.param = list(pt.cex = 0.5))
for (i in 1:n_spec) {
  segments(total_shape_within[all_links[,1],1,i], 
           total_shape_within[all_links[,1],2,i], 
           total_shape_within[all_links[,2],1,i], 
           total_shape_within[all_links[,2],2,i])
}

## ----outline landmarks--------------------------------------------------------
outline_lm <- papionin$outline$subset  # subset
k_out <- length(outline_lm)  # number of landmarks in the subset
outline_ar <- papionin_ar[outline_lm, , ]  # select the subset of landmark coordinates

## ----plot outline landmarks---------------------------------------------------
plot(papionin_ar[outline_lm,,1], pch = 16, col = "blue", asp = 1)
points(papionin_ar[-outline_lm,,1], pch = 4, col = "black", cex = 0.7)
text(papionin_ar[outline_lm,,1], label = outline_lm, pos = 1, col = "blue", cex = 0.8)
text(papionin_ar[-outline_lm,,1], label = c(1:k_out)[-outline_lm], pos = 1, col = "black", cex = 0.6)

## ----outline semilandmarks----------------------------------------------------
# Semilandmarks
outline_semi_lm <- papionin$outline$semi_lm
# Curves
outline_curves <- papionin$outline$curves
# Links between landmarks (for visualization)
outline_links <- papionin$outline$links

## ----plot outline fixed vs sliding landmarks----------------------------------
plot(outline_ar[,,1], pch = 16, col = "black", cex = 0.7, asp = 1)
points(outline_ar[-outline_semi_lm,,1], pch = 16, col = "red")
text(outline_ar[outline_semi_lm,,1], label = outline_semi_lm, pos = 1, col = "black", cex = 0.6)
text(outline_ar[-outline_semi_lm,,1], label = c(1:k_out)[-outline_semi_lm], pos = 1, col = "red", cex = 0.8)

## ----outline gpa sliding, message = FALSE-------------------------------------
outline_gpa <- Morpho::procSym(outline_ar, SMvector = outline_semi_lm, outlines = outline_curves, pcAlign = FALSE)  # sliding + GPA
outline_ar_slid <- outline_gpa$dataslide  # slid landmarks in the original space
outline_shape <- outline_gpa$rotated  # Procrustes coordinates
outline_shape_average <- outline_gpa$mshape  # average shape

## ----plot outline fixed vs slid landmarks-------------------------------------
plot(outline_ar_slid[,,1], pch = 16, col = "black", cex = 0.7, asp = 1)
points(outline_ar_slid[-outline_semi_lm,,1], pch = 16, col = "red")
text(outline_ar_slid[outline_semi_lm,,1], label = outline_semi_lm, pos = 1, col = "black", cex = 0.6)
text(outline_ar_slid[-outline_semi_lm,,1], label = c(1:k_out)[-outline_semi_lm], pos = 1, col = "red", cex = 0.8)

## ----plot outline shape-------------------------------------------------------
plotAllSpecimens(outline_shape, mean = FALSE, plot.param = list(pt.cex = 0.5))
for (i in 1:n_spec) {
  segments(outline_shape[outline_links[,1],1,i], 
           outline_shape[outline_links[,1],2,i], 
           outline_shape[outline_links[,2],1,i], 
           outline_shape[outline_links[,2],2,i])
}

## ----outline shape between species--------------------------------------------
# Computation of the average outline shape for each species
outline_shape_between <- array(NA, dim = c(k_out, 2, n_species))
for (i in 1:n_species) {
  spi <- which(species == levels(species)[i])
  outline_shape_between[,,i] <- mshape(outline_shape[, , spi])
}
# Plot outline shape between species
plotAllSpecimens(outline_shape_between, mean = FALSE, plot.param = list(pt.cex = 0.5))
for (i in 1:n_species) {
  segments(outline_shape_between[outline_links[,1],1,i], 
           outline_shape_between[outline_links[,1],2,i], 
           outline_shape_between[outline_links[,2],1,i], 
           outline_shape_between[outline_links[,2],2,i])
}

## ----outline shape within species---------------------------------------------
# Computation of the pooled individual within-species outline shape
outline_shape_within <- array(NA, dim = c(k_out, 2, n_spec))
for (i in 1:n_species) {
  spi <- which(species == levels(species)[i])
  for (j in 1:length(spi)) {
    outline_shape_within[,,spi[j]] <- outline_shape[, , spi[j]] - outline_shape_between[,,i] + mshape(outline_shape_between)
  }
}
# Plot the pooled individual within-species total shape
plotAllSpecimens(outline_shape_within, mean = FALSE, plot.param = list(pt.cex = 0.5))
for (i in 1:n_spec) {
  segments(outline_shape_within[outline_links[,1],1,i], 
           outline_shape_within[outline_links[,1],2,i], 
           outline_shape_within[outline_links[,2],1,i], 
           outline_shape_within[outline_links[,2],2,i])
}

## ----outline shape pca--------------------------------------------------------
# PCA
outline_shape_pca <- prcomp(two.d.array(outline_shape_between))
# Plot PC1 and PC2
col.sp <- c(rep("red", 3), rep("black", 2), rep("green", 2), rep ("blue", 8), "red", rep("green", 2))  # color vector
pch.sp <- c(rep(16, 4), 1, rep(16, 9), rep(1, 4))  # symbol vector
plot(outline_shape_pca$x[,1:2], las = 1, pch = pch.sp, col = col.sp, asp = 1, main = "Outline shape")
text(outline_shape_pca$x[,1:2], labels = levels(species), pos = 1, cex = 0.8)

## ----warping to the average outline-------------------------------------------
residual_shape <- tps.all(papionin_ar_slid, outline_ar_slid, outline_shape_average)  # warping to the average outline shape

## ----plot residual shape------------------------------------------------------
plotAllSpecimens(residual_shape, mean = FALSE, plot.param = list(pt.cex = 0.5))
for (i in 1:n_spec) {
  segments(residual_shape[all_links[,1],1,i], 
           residual_shape[all_links[,1],2,i], 
           residual_shape[all_links[,2],1,i], 
           residual_shape[all_links[,2],2,i])
}

## ----residual shape between species-------------------------------------------
# Computation of the average residual shape for each species
residual_shape_between <- array(NA, dim = c(k, 2, n_species))
for (i in 1:n_species) {
  spi <- which(species == levels(species)[i])
  residual_shape_between[,,i] <- mshape(residual_shape[, ,spi])
}
# Plot residual shape between species
plotAllSpecimens(residual_shape_between, mean = FALSE, plot.param = list(pt.cex = 0.5))
for (i in 1:n_species) {
  segments(residual_shape_between[all_links[,1],1,i], 
           residual_shape_between[all_links[,1],2,i], 
           residual_shape_between[all_links[,2],1,i], 
           residual_shape_between[all_links[,2],2,i])
}

## ----residual shape within species--------------------------------------------
# Computation of the pooled individual within-species residual shape
residual_shape_within <- array(NA, dim = c(k, 2, n_spec))
for (i in 1:n_species) {
  spi <- which(species == levels(species)[i])
  for (j in 1:length(spi)) {
    residual_shape_within[,,spi[j]] <- residual_shape[, , spi[j]] - residual_shape_between[,,i] + mshape(residual_shape_between)
  }
}
# Plot the pooled individual within-species residual shape
plotAllSpecimens(residual_shape_within, mean = FALSE, plot.param = list(pt.cex = 0.5))
for (i in 1:n_spec) {
  segments(residual_shape_within[all_links[,1],1,i], 
           residual_shape_within[all_links[,1],2,i], 
           residual_shape_within[all_links[,2],1,i], 
           residual_shape_within[all_links[,2],2,i])
}

## ----residual shape pca-------------------------------------------------------
# PCA
residual_shape_pca <- prcomp(two.d.array(residual_shape_between))
# Plot PC1 and PC2
plot(residual_shape_pca$x[,1:2], las = 1, pch = pch.sp, col = col.sp, asp = 1, main = "Residual shape")
text(residual_shape_pca$x[,1:2], labels = levels(species), pos = 1, cex = 0.8)

## ----partial warp decomposition-----------------------------------------------
d <- 50  # threshold for small-scale vs. large-scale components
papionin_be_pw <- create.pw.be(total_shape, total_shape_average, d)

## ----small-scale component----------------------------------------------------
small_shape_all <- xxyy.to.array(papionin_be_pw$Xsmall, p = k, k = 2) 
plot(small_shape_all[,1,], small_shape_all[,2,], asp = 1, las = 1, cex = 0.5, pch = 16, main = "Small-scale shape component", xlab = "X", ylab = "Y")
for (i in 1:n_spec) {
  segments(small_shape_all[all_links[,1],1,i], 
           small_shape_all[all_links[,1],2,i], 
           small_shape_all[all_links[,2],1,i], 
           small_shape_all[all_links[,2],2,i])
}

## ----large-scale component----------------------------------------------------
large_shape_all <- xxyy.to.array(papionin_be_pw$Xlarge, p = k, k = 2) 
plot(large_shape_all[,1,], large_shape_all[,2,], asp = 1, las = 1, cex = 0.5, pch = 16, main = "Large-scale shape component", xlab = "X", ylab = "Y")
for (i in 1:n_spec) {
  segments(large_shape_all[all_links[,1],1,i], 
           large_shape_all[all_links[,1],2,i], 
           large_shape_all[all_links[,2],1,i], 
           large_shape_all[all_links[,2],2,i])
}

## ----PW variance vs BE--------------------------------------------------------
# Computation of log BE^-1 for the (k-3) partial warps
logInvBE <- log((papionin_be_pw$bendingEnergy)^(-1))
# Computation of log PW variance for the (k-3) partial warps
logPWvar <- log(papionin_be_pw$variancePW)
# Linear regression of the log PW variance on the log BE^-1
mod <- lm(logPWvar ~ logInvBE)
# Plot log PW variance on log BE^-1 with regression line
plot(logInvBE, logPWvar, col = "white", asp = 1, main = "PW variance against inverse BE", xlab = "log 1/BE", ylab = "log PW variance", sub = paste("slope =", round(mod$coefficients[2], 2)))
text(logInvBE, logPWvar, labels = 1:(k-3), cex = 0.5)
abline(mod, col = "blue")

## ----partial warp decomposition between species-------------------------------
d <- 50  # threshold for small-scale vs. large-scale components
papionin_be_pw_between <- create.pw.be(total_shape_between, mshape(total_shape_between), d)

## ----small-scale component between species------------------------------------
small_shape_between <- xxyy.to.array(papionin_be_pw_between$Xsmall, p = k, k = 2) 
plot(small_shape_between[,1,], small_shape_between[,2,], asp = 1, las = 1, cex = 0.5, pch = 16, main = "Small-scale shape component", xlab = "X", ylab = "Y")
for (i in 1:n_species) {
  segments(small_shape_between[all_links[,1],1,i], 
           small_shape_between[all_links[,1],2,i], 
           small_shape_between[all_links[,2],1,i], 
           small_shape_between[all_links[,2],2,i])
}

## ----large-scale component between species------------------------------------
large_shape_between <- xxyy.to.array(papionin_be_pw_between$Xlarge, p = k, k = 2) 
plot(large_shape_between[,1,], large_shape_between[,2,], asp = 1, las = 1, cex = 0.5, pch = 16, main = "Large-scale shape component", xlab = "X", ylab = "Y")
for (i in 1:n_species) {
  segments(large_shape_between[all_links[,1],1,i], 
           large_shape_between[all_links[,1],2,i], 
           large_shape_between[all_links[,2],1,i], 
           large_shape_between[all_links[,2],2,i])
}

## ----PW variance vs BE between species----------------------------------------
# Computation of log BE^-1 for the (k-3) partial warps
logInvBE_between <- log((papionin_be_pw_between$bendingEnergy)^(-1))
# Computation of log PW variance for the (k-3) partial warps
logPWvar_between <- log(papionin_be_pw_between$variancePW)
# Linear regression of the log PW variance on the log BE^-1
mod_between <- lm(logPWvar_between ~ logInvBE_between)
# Plot log PW variance on log BE^-1 with regression line
plot(logInvBE_between, logPWvar_between, col = "white", asp = 1, main = "PW variance against inverse BE", xlab = "log 1/BE", ylab = "log PW variance", sub = paste("slope =", round(mod_between$coefficients[2], 2)))
text(logInvBE_between, logPWvar_between, labels = 1:(k-3), cex = 0.5)
abline(mod_between, col = "blue")

## ----partial warp decomposition within species--------------------------------
d <- 50  # threshold for small-scale vs. large-scale components
papionin_be_pw_within <- create.pw.be(total_shape_within, mshape(total_shape_within), d)

## ----small-scale component within species-------------------------------------
small_shape_within <- xxyy.to.array(papionin_be_pw_within$Xsmall, p = k, k = 2) 
plot(small_shape_within[,1,], small_shape_within[,2,], asp = 1, las = 1, cex = 0.5, pch = 16, main = "Small-scale shape component", xlab = "X", ylab = "Y")
for (i in 1:n_spec) {
  segments(small_shape_within[all_links[,1],1,i], 
           small_shape_within[all_links[,1],2,i], 
           small_shape_within[all_links[,2],1,i], 
           small_shape_within[all_links[,2],2,i])
}

## ----large-scale component within species-------------------------------------
large_shape_within <- xxyy.to.array(papionin_be_pw_within$Xlarge, p = k, k = 2) 
plot(large_shape_within[,1,], large_shape_within[,2,], asp = 1, las = 1, cex = 0.5, pch = 16, main = "Large-scale shape component", xlab = "X", ylab = "Y")
for (i in 1:n_spec) {
  segments(large_shape_within[all_links[,1],1,i], 
           large_shape_within[all_links[,1],2,i], 
           large_shape_within[all_links[,2],1,i], 
           large_shape_within[all_links[,2],2,i])
}

## ----PW variance vs BE within species-----------------------------------------
# Computation of log BE^-1 for the (k-3) partial warps
logInvBE_within <- log((papionin_be_pw_within$bendingEnergy)^(-1))
# Computation of log PW variance for the (k-3) partial warps
logPWvar_within <- log(papionin_be_pw_within$variancePW)
# Linear regression of the log PW variance on the log BE^-1
mod_within <- lm(logPWvar_within ~ logInvBE_within)
# Plot log PW variance on log BE^-1 with regression line
plot(logInvBE_within, logPWvar_within, col = "white", asp = 1, main = "PW variance against inverse BE", xlab = "log 1/BE", ylab = "log PW variance", sub = paste("slope =", round(mod_within$coefficients[2], 2)))
text(logInvBE_within, logPWvar_within, labels = 1:(k-3), cex = 0.5)
abline(mod_within, col = "blue")

