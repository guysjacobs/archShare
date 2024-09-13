### This script plots the simulated and observed summary statistics in summary statistic or PLS space, to assess overlap.
### A validation step to check that simulated summary statistics sufficiently capture that of the observed data.

# Import libraries
library(ggplot2)
library(gridExtra)

# Import data
setwd("/Users/hl636/Documents/Hirzi/Cambridge/Archaic sharing/ABC_archaic_working/")

for (j in seq(1,5)) {
  if (j==1) {
    pops <- "Nean_papua_noFrancoisUVBaining_easia"
    ABC_rej <- read.delim("Simulation_files/New_derivedStats_popMergeStats/ABCgrid_PapuaVsEAsia.T0=4058-T1=777.tightPriors100k.exponAncIntro.prod.allCorr.nean.power0.39.txt")
    ABC_obs<-read.delim("Observation_files/New_derivedStats_popMergeStats/Nean_papua_noFrancoisUVBaining_easia_wPopMergeStats.obs")
  } else if (j==2) {
    pops <- "Nean_papua_noFrancoisUVBaining_sasia"
    ABC_rej <- read.delim("Simulation_files/New_derivedStats_popMergeStats/ABCgrid_PapuaVsSAsia.T0=4054-T1=680.tightPriors100k.exponAncIntro.prod.allCorr.nean.power0.39.txt")
    ABC_obs<-read.delim("Observation_files/New_derivedStats_popMergeStats/Nean_papua_noFrancoisUVBaining_sasia_wPopMergeStats.obs")
  } else if (j==3) {
    pops <- "Nean_europe_easia"
    ABC_rej <- read.delim("Simulation_files/New_derivedStats_popMergeStats/ABCgrid_EuropeVsEAsia.T0=3792-T1=922.tightPriors100k.exponAncIntro.prod.allCorr.nean.power0.39.txt")
    ABC_obs<-read.delim("Observation_files/New_derivedStats_popMergeStats/Nean_europe_easia_wPopMergeStats.obs")
  } else if (j==4) {
    pops <- "Deni_papua_noFrancoisUVBaining_easia"
    ABC_rej <- read.delim("Simulation_files/New_derivedStats_popMergeStats/ABCgrid_PapuaVsEAsia.T0=4058-T1=777.tightPriors100k.exponNonPapIntro.prod.allCorr.deni.power0.29.txt")
    ABC_obs<-read.delim("Observation_files/New_derivedStats_popMergeStats/Deni_papua_noFrancoisUVBaining_easia_wPopMergeStats.obs")
  } else if (j==5) {
    pops <- "Deni_papua_noFrancoisUVBaining_sasia"
    ABC_rej <- read.delim("Simulation_files/New_derivedStats_popMergeStats/ABCgrid_PapuaVsSAsia.T0=4054-T1=680.tightPriors100k.exponNonPapIntro.prod.allCorr.deni.power0.29.txt")
    ABC_obs<-read.delim("Observation_files/New_derivedStats_popMergeStats/Deni_papua_noFrancoisUVBaining_sasia_wPopMergeStats.obs")
  }
  
  # Print population-pair names
  print(pops)

  # Select column with summary statistics / PLS LinearCombination_n
  ABC_rej<-ABC_rej[,c(10:ncol(ABC_rej))] # allCorr
  ABC_obs<-ABC_obs[,c(1:ncol(ABC_obs))]
  
  # Number of summary statistics or PLS components to plot
  num_PLS <- 24
  
  # Ensures an even number of summary statistics/PLS components. Since this plotter plots 2D plots, it requires an even number of summary statistics to plot. So in case of an odd number, just plot n+1 or n-1 summary statistics.
  if (num_PLS %% 2 !=0) {num_PLS <- num_PLS - 1}
  #if (num_PLS %% 2 !=0) {num_PLS <- num_PLS + 1}
  
  # For storing ggplot objects in a list using a loop: https://stackoverflow.com/questions/31993704/storing-ggplot-objects-in-a-list-from-within-loop-in-r
  plot_list <- list()
  #for (i in seq(1, (num_PLS/2 + 1))) {
  for (i in seq(1, (num_PLS/2))) {
    local({
      i <- i
      plot_list[[i]] <<- ggplot(data = ABC_rej, aes(x=ABC_rej[,2*i-1], y=ABC_rej[,2*i]) ) +
        geom_hex(bins = 35) + scale_fill_gradientn(colours=c("gray85","gray15"),name = "sim count",na.value=NA) +
        geom_hex(data = ABC_obs, bins = 70, aes(x=ABC_obs[,2*i-1], y=ABC_obs[,2*i], alpha=..count..), fill="red") +
        theme_bw() + xlab(colnames(ABC_obs)[2*i-1]) + ylab(colnames(ABC_obs)[2*i]) + labs("asd")
    })
  }
  
  grid.arrange(grobs = plot_list, ncol=3, top = paste0("Overlap of simulated & observed summary statistics - ", pops))
}
