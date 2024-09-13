# This script plots the marginal posterior densities (rejection sampling and GLM-regression) for the parameters of the archaic introgression-demographic model.

# Import libraries
library(MASS)
library(RColorBrewer)

# Define input variables
# Number of free parameters
nparams <- 6
# Number of observations; here 0-indexed so no of loci - 1
num_loci <- 0
# Plot Neanderthal or Denisovan ("nean" or "deni")
hom <- "deni"

# List of population-pair datasets
if (hom == "deni") {
  dataset_list <- c("Deni_papua_noFrancoisUVBaining_easia_wPopMergeStats", "Deni_papua_noFrancoisUVBaining_sasia_wPopMergeStats")
} else if (hom == "nean") {
  dataset_list <- c("Nean_papua_noFrancoisUVBaining_easia_wPopMergeStats", "Nean_papua_noFrancoisUVBaining_sasia_wPopMergeStats", "Nean_europe_easia_wPopMergeStats")
}

# Loop over population pairs
for (pop_pair in dataset_list) {
  
  # Print population pair
  print(pop_pair)
  
  # Define directory, prefix and suffix
  setwd(paste0("/Users/hl636/Documents/Hirzi/Cambridge/Archaic sharing/ABC_archaic_working/Estimation_results/New_derivedStats_popMergeStats/Estimation_results_numRetSims1000_diracPeakWidth0.005_jointPosteriors_prodRuns_allCor/Estimation_results_", pop_pair, "/"))
  prefix <- paste0("ABC_est_", pop_pair, "_model0_")
  suffix_marginals <- "MarginalPosteriorDensities_Obs"
  suffix_bestSims <- "BestSimsParamStats_Obs"
  
  # Plot the marginal posterior densities (refection sampling and GLM-regression) for all parameters
  custom_xlim_min <- c(0, 1000, 0, 1000, 0, 0)
  if (hom == "deni") {
    custom_xlim <- c(0.003, 1250, 0.055, 3250, 0.005, 1) # for Deni_prodRuns_allCor
    custom_ylim <- c(1500, 0.03, 155, 0.003, 850, 10) # for Deni_prodRuns_allCor
  } else if (hom == "nean") {
    custom_xlim <- c(0.0065, 3000, 0.0425, 3000, 0.0425, 1) # for Nean_prodRuns_allCor
    custom_ylim <- c(625, 0.007, 170, 0.008, 140, 2.5) # for Nean_prodRuns_allCor
  }

  par(mfrow=c(2,3))
  for(p in 1:nparams){
    # Plot all windows
    for (i in 0:num_loci) {
      # Should normalise rejection (and GLM) plot(s) first!
      # Read in data
      ABC_rej<-read.delim(paste(prefix, suffix_bestSims, i, ".txt", sep = ""))
      ABC_GLM<-read.delim(paste(prefix, suffix_marginals, i, ".txt", sep = ""))
      # Find modes of posterior distributions
      ABC_rej_density <- density(ABC_rej[,p+2])
      ABC_rej_density_df <- data.frame(ABC_rej_density$x, ABC_rej_density$y)
      colnames(ABC_rej_density_df) <- c("parameter_value", "density")
      mode_ABC_rej <- round(ABC_rej_density_df[,1][match(max(ABC_rej_density_df[,2]), ABC_rej_density_df[,2])], 5)
      mode_ABC_GLM <- round(ABC_GLM[,2*p][match(max(ABC_GLM[,2*p + 1]), ABC_GLM[,2*p + 1])], 5)
      # Plot
      #plot_label <- (paste(colnames(ABC_GLM)[2*p], mode_ABC_GLM, sep = ": "))
      plot_label <- (paste0(colnames(ABC_GLM)[2*p], ": ", mode_ABC_GLM, "(GLM), ", mode_ABC_rej, "(REJ)"))
      #plot(density(ABC_rej[,p+2]),main=pop_pair, type='l', lty=1, col='dodgerblue', ylim = c(0,max(ABC_GLM[[2*p + 1]])), lwd = 2, xlab = plot_label, ylab = NA, cex.lab = 1, cex.axis = 1, xaxs='i', yaxs='i')
      plot(density(ABC_rej[,p+2]),main=pop_pair, type='l', lty=1, col='dodgerblue', ylim = c(0,custom_ylim[p]), xlim = c(custom_xlim_min[p],custom_xlim[p]), lwd = 2, xlab = plot_label, ylab = NA, cex.lab = 1, cex.axis = 1, xaxs='i', yaxs='i')
      lines(ABC_GLM[[2*p]], ABC_GLM[[2*p + 1]], type='l',col='darkred', lwd = 2)
      abline(v=mode_ABC_rej, col="dodgerblue", lty=2)
      abline(v=mode_ABC_GLM, col="darkred", lty=2)
    }
  }
}
