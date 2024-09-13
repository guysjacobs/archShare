#!/usr/bin/env Rscript

### This script calculates and plot the residuals between spatially interpolated (raster) maps. 
### Here, the residuals reflect the frequency of archaic tag-SNPs identified in different regions compared to random SNPs matched for frequency in those regions.
### It uses the output of rspatial_kriging_vX.R as its input.

# Ran and tested with the following package versions 
"
> sessionInfo()
R version 4.2.1 (2022-06-23)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Monterey 12.7.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets 
[6] methods   base     

other attached packages:
[1] rnaturalearthdata_0.1.0 rnaturalearth_1.0.1    
[3] sf_1.0-8                stringr_1.5.0          
[5] leaflet_2.1.1           raster_3.6-20          
[7] sp_2.0-0                latticeExtra_0.6-30    
[9] RColorBrewer_1.1-3      gridExtra_2.3          
[11] rasterVis_0.51.2        lattice_0.20-45        

loaded via a namespace (and not attached):
[1] Rcpp_1.0.11        pillar_1.9.0      
[3] compiler_4.2.1     class_7.3-20      
[5] tools_4.2.1        digest_0.6.33     
[7] gtable_0.3.3       jsonlite_1.8.7    
[9] lifecycle_1.0.3    tibble_3.2.1      
[11] viridisLite_0.4.2  pkgconfig_2.0.3   
[13] png_0.1-8          rlang_1.1.1       
[15] DBI_1.1.3          cli_3.6.1         
[17] rstudioapi_0.14    crosstalk_1.2.0   
[19] parallel_4.2.1     hexbin_1.28.2     
[21] fastmap_1.1.1      interp_1.1-3      
[23] e1071_1.7-12       terra_1.7-39      
[25] dplyr_1.1.2        httr_1.4.6        
[27] htmlwidgets_1.6.2  generics_0.1.3    
[29] vctrs_0.6.3        classInt_0.4-8    
[31] grid_4.2.1         tidyselect_1.2.0  
[33] glue_1.6.2         R6_2.5.1          
[35] jpeg_0.1-10        fansi_1.0.4       
[37] deldir_1.0-6       magrittr_2.0.3    
[39] htmltools_0.5.5    codetools_0.2-18  
[41] units_0.8-0        KernSmooth_2.23-20
[43] utf8_1.2.3         stringi_1.7.12    
[45] proxy_0.4-27       zoo_1.8-11 
"

# Load libraries
library("rasterVis")
library("gridExtra")
library("RColorBrewer")
library("lattice")
library("latticeExtra")
library("raster")
library("leaflet")
library("stringr")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

# Set working directory and run parameters
setwd("/Users/hl636/Documents/Hirzi/Cambridge/Archaic sharing/")
output_rasters <-stack("Archaic_maps_raw_revised_2024")
pops_to_compare <- c("deniHCSS35unique.subset_papuaExclFrancoisUVBaining.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc.w5000.noAfPoly", "deniHCSS35unique.subset_papuaExclFrancoisUVBaining.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc-resample0.w5000.noAfPoly", "deniHCSS35unique.continent_easia.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc.w5000.noAfPoly", "deniHCSS35unique.continent_easia.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc-resample0.w5000.noAfPoly", "deniHCSS35unique.continent_sasia.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc.w5000.noAfPoly", "deniHCSS35unique.continent_sasia.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc-resample0.w5000.noAfPoly", "neanHCSS35unique.subset_papuaExclFrancoisUVBaining.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc.w5000.noAfPoly", "neanHCSS35unique.subset_papuaExclFrancoisUVBaining.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc-resample0.w5000.noAfPoly", "neanHCSS35unique.continent_easia.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc.w5000.noAfPoly", "neanHCSS35unique.continent_easia.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc-resample0.w5000.noAfPoly", "neanHCSS35unique.continent_sasia.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc.w5000.noAfPoly", "neanHCSS35unique.continent_sasia.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc-resample0.w5000.noAfPoly", "neanHCSS35unique.continent_europe.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc.w5000.noAfPoly", "neanHCSS35unique.continent_europe.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc-resample0.w5000.noAfPoly") # final ascertainment scheme
pops_to_compare_names <- c("Deni_Papua", "Deni_Papua_Resample", "Deni_EAsia", "Deni_EAsia_Resample", "Deni_SAsia", "Deni_SAsia_Resample", "Nean_Papua", "Nean_Papua_Resample", "Nean_EAsia", "Nean_EAsia_Resample", "Nean_SAsia", "Nean_SAsia_Resample", "Nean_Europe", "Nean_Europe_Resample")
# Set number of colour breaks for plots
nColBreaks <- 17 # this will generate nColBreaks-1 colours
# Set final resolution of raster plots
res_final <- 1e6

# Add points and labels
cols <- c(rep("gray30", 12), "black") # consistent colour palette with other figures
cols2 <- c(rep("black", 12), "gray50")
cex_vec <- c(rep(0.85,12), 0.35)
pch_vec <- c(rep(21,12), 19)
coords_wContinentInfo = read.table("/Users/hl636/Documents/Hirzi/Cambridge/Archaic sharing/Individual_GPS_locations_withContinent_revised.txt", col.names = c("IND", "LATITUDE", "LONGITUDE","CONTINENT"), check.names=FALSE)
coords_wContinentInfo$col <- as.numeric(factor(coords_wContinentInfo$CONTINENT))
inds_to_remove_idx <- which(str_extract(coords_wContinentInfo$IND, "^.{2}") == "UV") # remove Baining samples
coords_wContinentInfo <- coords_wContinentInfo[-inds_to_remove_idx, ]

# Define coastline
coast = rnaturalearth::ne_coastline(returnclass = "sv")
crs(coast) <- "+proj=longlat +datum=WGS84"
coast <- as(coast, "Spatial")

#Plot residuals (substraction)
pos_residuals_stack <- stack()
neg_residuals_stack <- stack()

for (i in seq(1, length(pops_to_compare)/2)) {
  idx <- (2*i)-1
  print(idx)
  pos_residuals_raster <- output_rasters[[idx]]-output_rasters[[idx+1]]
  names(pos_residuals_raster) <- pops_to_compare_names[idx]
  neg_residuals_raster <- output_rasters[[idx+1]]-output_rasters[[idx]]
  names(neg_residuals_raster) <- pops_to_compare_names[idx]
  pos_residuals_stack <- stack(pos_residuals_stack, pos_residuals_raster)
  neg_residuals_stack <- stack(neg_residuals_stack, neg_residuals_raster)
}

# Make empty list
levelplot_list <- list()
levelplot_fixedz_list <- list()
coords_wContinentInfo_pop <- coords_wContinentInfo
for (j in 1:nlayers(pos_residuals_stack)) {
  # Define coordinate colours
  coords_wContinentInfo_pop[, ncol(coords_wContinentInfo_pop)+1] <- coords_wContinentInfo_pop$col
  colnames(coords_wContinentInfo_pop)[ncol(coords_wContinentInfo_pop)] <- paste0("col", j)
  if(str_split(names(pos_residuals_stack)[j], "_")[[1]][2] == "Papua") {
    coords_wContinentInfo_pop[, ncol(coords_wContinentInfo_pop)][coords_wContinentInfo_pop$CONTINENT != "papua"] <- 13
  } else if(str_split(names(pos_residuals_stack)[j], "_")[[1]][2] == "EAsia") {
    coords_wContinentInfo_pop[, ncol(coords_wContinentInfo_pop)][coords_wContinentInfo_pop$CONTINENT != "easia"] <- 13
  } else if(str_split(names(pos_residuals_stack)[j], "_")[[1]][2] == "SAsia") {
    coords_wContinentInfo_pop[, ncol(coords_wContinentInfo_pop)][coords_wContinentInfo_pop$CONTINENT != "sasia"] <- 13
  } else if(str_split(names(pos_residuals_stack)[j], "_")[[1]][2] == "Europe") {
    coords_wContinentInfo_pop[, ncol(coords_wContinentInfo_pop)][coords_wContinentInfo_pop$CONTINENT != "europe"] <- 13
  }
  fd2 <- coords_wContinentInfo_pop
  coordinates(fd2) <- ~LONGITUDE + LATITUDE
  
  # Plot residuals on diverging colour gradient centred on zero.
  intr_minus_rand <- pos_residuals_stack[[j]]
  rand_minus_intr <- neg_residuals_stack[[j]]
  intr_minus_rand[intr_minus_rand < 0] <- NA
  rand_minus_intr[intr_minus_rand < 0] <- NA
  rand_minus_intr.masked <- mask(rand_minus_intr, intr_minus_rand, inverse=TRUE)
  intr_minus_rand.masked <- mask(intr_minus_rand, rand_minus_intr)
  # Make one raster all positive values and the other all negative values, so that we can combine the two into one plot (raster)
  rand_minus_intr.masked.negative <- 0 - rand_minus_intr.masked
  # Use the 'cover' function to overlay the two rasters together
  current.continuous.turnover <- cover(rand_minus_intr.masked.negative, intr_minus_rand.masked)
  # Project back to longlat WGS84 projection
  current.continuous.turnover <- projectRaster(current.continuous.turnover, crs='+proj=longlat +datum=WGS84')
  e <- c(-27.5, 178.5148, -48.17816, 78.04584)
  current.continuous.turnover <- crop(current.continuous.turnover, e)
  # Output free z-axis plots
  if (j==1) {
    current.continuous.turnover.plot <- levelplot(current.continuous.turnover,par.settings=RdBuTheme(), at=seq(-max(abs(cellStats(current.continuous.turnover, range))), max(abs(cellStats(current.continuous.turnover, range))), len=nColBreaks),  maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col1], pch = pch_vec[coords_wContinentInfo_pop$col1], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col1], col = cols2[coords_wContinentInfo_pop$col1])) +
      layer(sp.lines(coast))
  } else if (j==2) {
    current.continuous.turnover.plot <- levelplot(current.continuous.turnover,par.settings=RdBuTheme(), at=seq(-max(abs(cellStats(current.continuous.turnover, range))), max(abs(cellStats(current.continuous.turnover, range))), len=nColBreaks),  maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col2], pch = pch_vec[coords_wContinentInfo_pop$col2], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col2], col = cols2[coords_wContinentInfo_pop$col2])) +
      layer(sp.lines(coast))
  } else if (j==3) {
    current.continuous.turnover.plot <- levelplot(current.continuous.turnover,par.settings=RdBuTheme(), at=seq(-max(abs(cellStats(current.continuous.turnover, range))), max(abs(cellStats(current.continuous.turnover, range))), len=nColBreaks),  maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col3], pch = pch_vec[coords_wContinentInfo_pop$col3], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col3], col = cols2[coords_wContinentInfo_pop$col3])) +
      layer(sp.lines(coast))
  } else if (j==4) {
    current.continuous.turnover.plot <- levelplot(current.continuous.turnover,par.settings=RdBuTheme(), at=seq(-max(abs(cellStats(current.continuous.turnover, range))), max(abs(cellStats(current.continuous.turnover, range))), len=nColBreaks),  maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col4], pch = pch_vec[coords_wContinentInfo_pop$col4], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col4], col = cols2[coords_wContinentInfo_pop$col4])) +
      layer(sp.lines(coast))
  } else if (j==5) {
    current.continuous.turnover.plot <- levelplot(current.continuous.turnover,par.settings=RdBuTheme(), at=seq(-max(abs(cellStats(current.continuous.turnover, range))), max(abs(cellStats(current.continuous.turnover, range))), len=nColBreaks),  maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col5], pch = pch_vec[coords_wContinentInfo_pop$col5], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col5], col = cols2[coords_wContinentInfo_pop$col5])) +
      layer(sp.lines(coast))
  } else if (j==6) {
    current.continuous.turnover.plot <- levelplot(current.continuous.turnover,par.settings=RdBuTheme(), at=seq(-max(abs(cellStats(current.continuous.turnover, range))), max(abs(cellStats(current.continuous.turnover, range))), len=nColBreaks),  maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col6], pch = pch_vec[coords_wContinentInfo_pop$col6], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col6], col = cols2[coords_wContinentInfo_pop$col6])) +
      layer(sp.lines(coast))
  } else if (j==7) {
    current.continuous.turnover.plot <- levelplot(current.continuous.turnover,par.settings=RdBuTheme(), at=seq(-max(abs(cellStats(current.continuous.turnover, range))), max(abs(cellStats(current.continuous.turnover, range))), len=nColBreaks),  maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col7], pch = pch_vec[coords_wContinentInfo_pop$col7], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col7], col = cols2[coords_wContinentInfo_pop$col7])) +
      layer(sp.lines(coast))
  }
  # Output fixed z-axis plots
  ## Make a asymmetric diverging colour gradient centred around 0. For reference, see: https://stackoverflow.com/questions/49126405/how-to-set-up-asymmetrical-color-gradient-for-a-numerical-variable-in-leaflet-in and https://stackoverflow.com/questions/29262824/r-center-color-palette-on-0
  print(pops_to_compare_names[2*j-1])
  print(cellStats(current.continuous.turnover, range))
  my.at.data <- seq(cellStats(current.continuous.turnover, range)[1], cellStats(current.continuous.turnover, range)[2], len=nColBreaks)
  my.at.nean <- seq(-0.0147, 0.0423, len=nColBreaks) # full range, includes some arctic islands
  my.at.deni <- seq(-0.02, 0.036, len=nColBreaks) # full range, includes new caledonia, vanuatu, east solomon islands
  if (j <= 3) {
    my.at <- my.at.deni
  } else if (j > 3) {
    my.at <- my.at.nean
  }
  myColorkey <- list(at=my.at)
  ## Make vector of colors for negative values
  rc1 <- colorRampPalette(colors = c("#B2182B", "#F4AA88", "#F7F7F7"), space = "Lab")(length(which(my.at<0)))
  ## Make vector of colors for positive values
  rc2 <- colorRampPalette(colors = c("#F7F7F7", "#98C8DF", "#2166AC"), space = "Lab")(length(which(my.at>0)))
  ## Combine the two color palettes
  rampcols <- c(rc1, rc2)
  myTheme <- rasterTheme(region = rampcols)
  # And plot
  if (j==1) {
    current.continuous.turnover.fixedz.plot <- levelplot(current.continuous.turnover, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col1], pch = pch_vec[coords_wContinentInfo_pop$col1], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col1], col = cols2[coords_wContinentInfo_pop$col1])) +
      layer(sp.lines(coast))
  } else if (j==2) {
    current.continuous.turnover.fixedz.plot <- levelplot(current.continuous.turnover, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col2], pch = pch_vec[coords_wContinentInfo_pop$col2], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col2], col = cols2[coords_wContinentInfo_pop$col2])) +
      layer(sp.lines(coast))
  } else if (j==3) {
    current.continuous.turnover.fixedz.plot <- levelplot(current.continuous.turnover, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col3], pch = pch_vec[coords_wContinentInfo_pop$col3], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col3], col = cols2[coords_wContinentInfo_pop$col3])) +
      layer(sp.lines(coast))
  } else if (j==4) {
    current.continuous.turnover.fixedz.plot <- levelplot(current.continuous.turnover, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col4], pch = pch_vec[coords_wContinentInfo_pop$col4], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col4], col = cols2[coords_wContinentInfo_pop$col4])) +
      layer(sp.lines(coast))
  } else if (j==5) {
    current.continuous.turnover.fixedz.plot <- levelplot(current.continuous.turnover, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col5], pch = pch_vec[coords_wContinentInfo_pop$col5], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col5], col = cols2[coords_wContinentInfo_pop$col5])) +
      layer(sp.lines(coast))
  } else if (j==6) {
    current.continuous.turnover.fixedz.plot <- levelplot(current.continuous.turnover, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col6], pch = pch_vec[coords_wContinentInfo_pop$col6], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col6], col = cols2[coords_wContinentInfo_pop$col6])) +
      layer(sp.lines(coast))
  } else if (j==7) {
    current.continuous.turnover.fixedz.plot <- levelplot(current.continuous.turnover, par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = pops_to_compare_names[2*j-1], legend=list(top=list(fun=grid::textGrob("Residuals", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col7], pch = pch_vec[coords_wContinentInfo_pop$col7], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col7], col = cols2[coords_wContinentInfo_pop$col7])) +
      layer(sp.lines(coast))
  }
  # Store in list
  levelplot_list[[pops_to_compare_names[2*j-1]]] <- current.continuous.turnover.plot
  levelplot_fixedz_list[[pops_to_compare_names[2*j-1]]] <- current.continuous.turnover.fixedz.plot
}

# Plot in grid. Export 20x12 (PDF), 2000x1200 (tiff)
grid.arrange(levelplot_list[[1]], levelplot_list[[2]], levelplot_list[[3]], ncol=2)
grid.arrange(levelplot_list[[4]], levelplot_list[[5]], levelplot_list[[6]], levelplot_list[[7]], ncol=2)

## Plot predictions
# Project back to longlat WGS84 projection
output_rasters <- projectRaster(output_rasters, crs='+proj=longlat +datum=WGS84')
e <- c(-27.5, 178.5148, -48.17816, 78.04584)
output_rasters <- crop(output_rasters, e)
# And plot
levelplot_predictions_list <- list()
for (i in seq(1, length(pops_to_compare)/2)) {
  idx <- (2*i)-1
  print(names(output_rasters)[idx])
  # Plotting parameters
  my.at <- seq(cellStats(output_rasters[[idx]], range)[1], cellStats(output_rasters[[idx]], range)[2], len=nColBreaks) # free z-axes
  colPal <- brewer.pal('YlOrRd', n = 9) #n=9 here has no bearing on number of colour breaks in plot as colours are interpolated via myTheme (via my.at)
  myColorkey <- list(at=my.at)
  myTheme <- rasterTheme(region = colPal)
  # And plot
  if (i==1) {
    prediction.plot <- levelplot(output_rasters[[idx]], par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = names(output_rasters)[idx], legend=list(top=list(fun=grid::textGrob("Introgression rate", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col1], pch = pch_vec[coords_wContinentInfo_pop$col1], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col1], col = cols2[coords_wContinentInfo_pop$col1])) +
      layer(sp.lines(coast))
  } else if (i==2) {
    prediction.plot <- levelplot(output_rasters[[idx]], par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = names(output_rasters)[idx], legend=list(top=list(fun=grid::textGrob("Introgression rate", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col2], pch = pch_vec[coords_wContinentInfo_pop$col2], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col2], col = cols2[coords_wContinentInfo_pop$col2])) +
      layer(sp.lines(coast))
  } else if (i==3) {
    prediction.plot <- levelplot(output_rasters[[idx]], par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = names(output_rasters)[idx], legend=list(top=list(fun=grid::textGrob("Introgression rate", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col3], pch = pch_vec[coords_wContinentInfo_pop$col3], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col3], col = cols2[coords_wContinentInfo_pop$col3])) +
      layer(sp.lines(coast))
  } else if (i==4) {
    prediction.plot <- levelplot(output_rasters[[idx]], par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = names(output_rasters)[idx], legend=list(top=list(fun=grid::textGrob("Introgression rate", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col4], pch = pch_vec[coords_wContinentInfo_pop$col4], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col4], col = cols2[coords_wContinentInfo_pop$col4])) +
      layer(sp.lines(coast))
  } else if (i==5) {
    prediction.plot <- levelplot(output_rasters[[idx]], par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = names(output_rasters)[idx], legend=list(top=list(fun=grid::textGrob("Introgression rate", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col5], pch = pch_vec[coords_wContinentInfo_pop$col5], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col5], col = cols2[coords_wContinentInfo_pop$col5])) +
      layer(sp.lines(coast))
  } else if (i==6) {
    prediction.plot <- levelplot(output_rasters[[idx]], par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = names(output_rasters)[idx], legend=list(top=list(fun=grid::textGrob("Introgression rate", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col6], pch = pch_vec[coords_wContinentInfo_pop$col6], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col6], col = cols2[coords_wContinentInfo_pop$col6])) +
      layer(sp.lines(coast))
  } else if (i==7) {
    prediction.plot <- levelplot(output_rasters[[idx]], par.settings = myTheme, at=my.at, colorkey=myColorkey, maxpixels = res_final, margin = FALSE, main = names(output_rasters)[idx], legend=list(top=list(fun=grid::textGrob("Introgression rate", y=1, x=1.06)))) + layer(sp.points(fd2, cex = cex_vec[coords_wContinentInfo_pop$col7], pch = pch_vec[coords_wContinentInfo_pop$col7], lwd = 0.5, fill = cols[coords_wContinentInfo_pop$col7], col = cols2[coords_wContinentInfo_pop$col7])) +
      layer(sp.lines(coast))
  }
  # Store in list
  levelplot_predictions_list[[names(output_rasters)[idx]]] <- prediction.plot
}

# Plot in grid. Export 20x12 (PDF), 2000x1200 (tiff)
grid.arrange(levelplot_predictions_list[[1]], levelplot_predictions_list[[2]], levelplot_predictions_list[[3]], ncol=2)
grid.arrange(levelplot_predictions_list[[4]], levelplot_predictions_list[[5]], levelplot_predictions_list[[6]], levelplot_predictions_list[[7]], ncol=2)
