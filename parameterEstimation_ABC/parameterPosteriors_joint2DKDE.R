# This script plots 2DKDE posteriors for joint (ly estimated) parameters

# Import libraries
library(MASS)
library(RColorBrewer)
library(MASS)

# Read in ABCEstimation data
setwd("/Users/hl636/Documents/Hirzi/Cambridge/Archaic sharing/ABC_archaic_working/Estimation_results/New_derivedStats_popMergeStats/Estimation_results_numRetSims1000_diracPeakWidth0.005_jointPosteriors_prodRuns_allCor/Estimation_results_Deni_papua_noFrancoisUVBaining_easia_wPopMergeStats/")
prefix <- "ABC_est_Deni_papua_noFrancoisUVBaining_easia_wPopMergeStats_model0_jointPosterior_3_4_Obs" #I_prop_a___I_Tb
#prefix <- "ABC_est_Deni_papua_noFrancoisUVBaining_easia_wPopMergeStats_model0_jointPosterior_1_3_Obs" #I_prop_x___I_prop_a
#prefix <- "ABC_est_Deni_papua_noFrancoisUVBaining_easia_wPopMergeStats_model0_jointPosterior_1_5_Obs" #I_prop_x___I_prop_b
ABC_GLM<-read.delim(paste(prefix, 0, ".txt", sep = ""))

# Acquire prior range for joint parameters
prior_jointParam1 <- c(round(min(ABC_GLM[2]), digits = 3),round(max(ABC_GLM[2]), digits = 3))
prior_jointParam2 <- c(round(min(ABC_GLM[3]), digits = 3),round(max(ABC_GLM[3]), digits = 3))
input_data <- ABC_GLM

# For the joint posteriors, acquire the size of the grid along one axis (e.g. the number of sampled points per parameter). Assuming square grid, we can simply take the square root of the number of samples.
density_points = sqrt(nrow(ABC_GLM))
# Let's calculate symmetry of posterior distribution (about the diagonal)
# Recall that there are points above the diagonal, below the diagonal and on the diagonal. To reduce degrees of freedom, we first account for the points on the diagonal (i.e. where mLH=Mhl).
symmetry <- sum((ABC_GLM[,2] == ABC_GLM[,3])*ABC_GLM$density)
# Then the integral of the area above the diagonal equals one minus the integral of the area below the diagonal. So we can just output one variable.
asymmetry <- sum((ABC_GLM[,2] < ABC_GLM[,3])*ABC_GLM$density) / (sum(ABC_GLM$density) - symmetry)
# Then prepare the density plot
dens <- list()
# This (x-axis) refers to second column of ABC_GLM
dens[[1]] <- seq(prior_jointParam1[1],prior_jointParam1[2], length.out = density_points)
# This (y-axis) refers to third column of ABC_GLM
dens[[2]] <- seq(prior_jointParam2[1],prior_jointParam2[2], length.out = density_points)
# We transform the density vector into matrix.
dens[[3]] <- matrix(ABC_GLM$density,nrow=density_points, ncol = density_points)
names(dens) <- c("x", "y", "z")

# Define argument behaviour
lim_x <- c(prior_jointParam1[1],prior_jointParam1[2]);
lim_y <- c(prior_jointParam2[1],prior_jointParam2[2]);

# Plot
zlim <- (max(dens$z) * 1.1)
par(mfrow=c(1,1),mar=c(5.1, 4.1, 4.1, 2.1))
filled.contour(dens, xlim=lim_x, ylim=lim_y, xlab=colnames(ABC_GLM)[2], ylab=colnames(ABC_GLM)[3], cex.lab = 1, levels=seq(0, zlim, length.out = 40), 
               color.palette=colorRampPalette(c('white','red3')))
