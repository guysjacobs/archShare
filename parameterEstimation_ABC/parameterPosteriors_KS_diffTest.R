# This scripts performs Kolmogorov-Smirnov Tests between marginal posteriors, i.e. to test if two marginal posterior (parameter) distributions are significantly different from each other.
# For reference on the K-S test, see: https://stats.stackexchange.com/questions/471732/intuitive-explanation-of-kolmogorov-smirnov-test

# Define population list
dataset_list <- c("Nean_papua_noFrancoisUVBaining_easia_wPopMergeStats", "Nean_papua_noFrancoisUVBaining_sasia_wPopMergeStats", "Nean_europe_easia_wPopMergeStats", "Deni_papua_noFrancoisUVBaining_easia_wPopMergeStats", "Deni_papua_noFrancoisUVBaining_sasia_wPopMergeStats")
pop_names <- c("Nean_papua_easia", "Nean_papua_sasia", "Nean_europe_easia", "Deni_papua_easia", "Deni_papua_sasia")

# Intialise counter (only for plot labels)
i <- 1
# Initialise  results vectors
ks_time_vec <- list()
ks_prop_vec <- list()

# Loop over population pairs
for (pop_pair in dataset_list) {
  
  # Print population pair
  print(pop_pair)

  # Define working directory and files
  setwd(paste0("/Users/hl636/Documents/Hirzi/Cambridge/Archaic sharing/ABC_archaic_working/Estimation_results/New_derivedStats_popMergeStats/Estimation_results_numRetSims1000_diracPeakWidth0.005_jointPosteriors_prodRuns_allCor/Estimation_results_", pop_pair, "/"))
  prefix <- paste0("ABC_est_", pop_pair, "_model0_")
  suffix_bestSims <- "BestSimsParamStats_Obs"

  # Read in top sims
  ABC_rej<-read.delim(paste(prefix, suffix_bestSims, 0, ".txt", sep = ""))

  # Perform K-S test. Alternative null hypothesis can be one of "less", "greater", "two.sided"
  ks_time <- ks.test(ABC_rej$I_Ta, ABC_rej$I_Tb, alternative = "two.sided", exact = TRUE)
  #ks_time <- ks.test(ABC_rej$I_Ta, ABC_rej$I_Tb, alternative = "l", exact = TRUE)
  ks_prop <- ks.test(ABC_rej$I_prop_a, ABC_rej$I_prop_b, alternative = "two.sided", exact = TRUE)
  #ks_prop <- ks.test(ABC_rej$I_prop_a, ABC_rej$I_prop_b, alternative = "l", exact = TRUE)
  #ks_prop <- ks.test(ABC_rej$I_prop_a, ABC_rej$I_prop_b, alternative = "l") 
  
  # Is p-value = NA because there is no overlap? Because with no overlap, the KS statistic which is the max y-value difference between the two ecdfs can't be calculated.
  # If so, then theoretically they are totally distinct, i.e. p-value converges to 0.
  
  # Append K-S results to vector
  ks_time_vec[[i]] <- ks_time
  ks_prop_vec[[i]] <- ks_prop
  
  # Plot ECDFs to check
  par(mfrow=c(1,2))
  plot(ecdf(ABC_rej$I_Ta), col='dodgerblue', xlim = range(c(ABC_rej$I_Ta, ABC_rej$I_Tb)), main=paste0(pop_names[i], "_iTa_vs_iTb"))
  plot(ecdf(ABC_rej$I_Tb), col='darkred', add = TRUE)
  plot(ecdf(ABC_rej$I_prop_a), col='dodgerblue', xlim = range(c(ABC_rej$I_prop_a, ABC_rej$I_prop_b)), main=paste0(pop_names[i], "_iPropA_vs_iPropB"))
  plot(ecdf(ABC_rej$I_prop_b), col='darkred', add = TRUE)
  
  # Update counter
  i <- i + 1
}
