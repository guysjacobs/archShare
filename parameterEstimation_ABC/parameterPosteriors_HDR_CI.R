##### This script estimates modes and credible intervals of parameter posteriors #####

## Define variables
# Number of free parameters
nparams <- 5
# Highest density credible interval. E.g. 95%
cred_interval <- 0.95
# Population pairs
dataset_list <- c("Nean_papua_noFrancoisUVBaining_easia_wPopMergeStats", "Nean_papua_noFrancoisUVBaining_sasia_wPopMergeStats", "Nean_europe_easia_wPopMergeStats", "Deni_papua_noFrancoisUVBaining_easia_wPopMergeStats", "Deni_papua_noFrancoisUVBaining_sasia_wPopMergeStats")

# Loop over population pairs
HDIs_summary <- list()
modes_summary <- list()
for (pop_pair in dataset_list) {
  ## Read in data
  # Change working directory accordingly
  setwd(paste0("/Users/hl636/Documents/Hirzi/Cambridge/Archaic sharing/ABC_archaic_working/Estimation_results/New_derivedStats_popMergeStats/Estimation_results_numRetSims1000_diracPeakWidth0.005_jointPosteriors_prodRuns_allCor/Estimation_results_", pop_pair, "/"))
  prefix <- paste0("ABC_est_", pop_pair, "_model0_")
  suffix_bestSims <- "BestSimsParamStats_Obs"
  ABC_rej<-read.delim(paste(prefix, suffix_bestSims, "0.txt", sep = ""))
  ## Loop over all parameter posteriors
  modes <- list()
  HDIs <- list()
  for(p in 1:nparams){
    # Parameter name
    param_name <- colnames(ABC_rej)[3:7][p]
    # Kernel density estimation
    ABC_rej_density <- density(ABC_rej[,p+2])
    ABC_rej_density_df <- data.frame(ABC_rej_density$x, ABC_rej_density$y)
    colnames(ABC_rej_density_df) <- c("parameter_value", "density")
    # Find mode
    mode_ABC_rej <- round(ABC_rej_density_df[,1][match(max(ABC_rej_density_df[,2]), ABC_rej_density_df[,2])], 5)
    modes[[param_name]] <- mode_ABC_rej
    # Calculate credible interval (here: highest density interval)
    const2d <- sum(ABC_rej_density_df$density)
    spxx2d <- sort(ABC_rej_density_df$density, decreasing = TRUE) / const2d
    HDI_dens <- spxx2d[which(cumsum(spxx2d) >= cred_interval[1])[1]] * const2d
    ABC_rej_density_df[,ncol(ABC_rej_density_df)+1] <- ABC_rej_density_df$density - HDI_dens
    colnames(ABC_rej_density_df)[ncol(ABC_rej_density_df)] <- "delta"
    ABC_rej_density_filtered <- ABC_rej_density_df[ABC_rej_density_df$delta >= 0,]
    #HDI_min <- min(ABC_rej_density_filtered$parameter_value)
    HDI_min <- max(0, min(ABC_rej_density_filtered$parameter_value)) # NOTE: kernel density  estimation can lead to negative values even if prior range is defined for positive values only. As negative values are undefined, we collapse them to 0.
    HDI_max <- max(ABC_rej_density_filtered$parameter_value)
    HDI <- c(HDI_min, HDI_max)
    HDIs[[param_name]] <- HDI
    # Plot
    plot(density(ABC_rej[,p+2]),main=pop_pair, type='l', lty=1, col='dodgerblue', lwd = 2, xlab = param_name, ylab = NA, cex.lab = 1, cex.axis = 1, xaxs='i', yaxs='i')
    abline(v=mode_ABC_rej, col="darkred", lty=2, lwd = 2)
    abline(v=HDI_min, col="dodgerblue", lty=2, lwd = 2)
    abline(v=HDI_max, col="dodgerblue", lty=2, lwd = 2)
  }
  HDIs_summary[[pop_pair]] <- HDIs
  modes_summary[[pop_pair]] <- modes
}

print(HDIs_summary)
print(modes_summary)
