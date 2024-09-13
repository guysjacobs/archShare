#!/usr/bin/env Rscript

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
[1] raster_3.6-20       sp_2.0-0           
[3] stringi_1.7.12      terra_1.7-39       
[5] sf_1.0-8            automap_1.1-9      
[7] gstat_2.0-9         rnaturalearth_1.0.1
[9] optparse_1.7.3     

loaded via a namespace (and not attached):
[1] Rcpp_1.0.11             plyr_1.8.8             
[3] pillar_1.9.0            compiler_4.2.1         
[5] class_7.3-20            tools_4.2.1            
[7] xts_0.13.1              gtable_0.3.3           
[9] jsonlite_1.8.7          lifecycle_1.0.3        
[11] tibble_3.2.1            lattice_0.20-45        
[13] pkgconfig_2.0.3         rlang_1.1.1            
[15] DBI_1.1.3               cli_3.6.1              
[17] rstudioapi_0.14         parallel_4.2.1         
[19] rnaturalearthdata_0.1.0 e1071_1.7-12           
[21] dplyr_1.1.2             httr_1.4.6             
[23] generics_0.1.3          vctrs_0.6.3            
[25] classInt_0.4-8          grid_4.2.1             
[27] tidyselect_1.2.0        getopt_1.20.3          
[29] reshape_0.8.9           spacetime_1.2-8        
[31] glue_1.6.2              R6_2.5.1               
[33] fansi_1.0.4             ggplot2_3.4.2          
[35] magrittr_2.0.3          scales_1.2.1           
[37] stars_0.5-6             intervals_0.15.2       
[39] codetools_0.2-18        units_0.8-0            
[41] abind_1.4-5             colorspace_2.1-0       
[43] utf8_1.2.3              KernSmooth_2.23-20     
[45] proxy_0.4-27            munsell_0.5.0          
[47] lwgeom_0.2-9            FNN_1.1.3.1            
[49] zoo_1.8-11
"

# Load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(rnaturalearth))
suppressPackageStartupMessages(library(gstat))
suppressPackageStartupMessages(library(automap))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(terra))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(raster))

### This program uses gstat (and automap) to do spatial kriging on e.g. frequency of other data.
### There are a few options for flexibility. It would be pretty easy and potentially interesting
### to adapt for e.g. co-kriging, or for other forms of spatial analysis.

option_list = list(
  make_option(c("-l", "--location_file"), type="character", default="/home/guy/Dropbox/Transfer/Indonesia_Diversity175_data/geography/Individual_GPS_locations.txt", 
              help="name of file containing locations. These must have columns IND,LATITUDE,LONGITUDE, or namese LATITUDE,LONGITUDE and be in the same order as the freq file.", metavar="character"),
  make_option(c("-f", "--freq_file"), type="character", default="/home/guy/Dropbox/Transfer/chunk_analyses/coverage_countsV4_manycombined_totbp_inds.csv", 
              help="name of file containing frequencies. These must either have individual IDs in a column Ind or be in the same order as the location file.", metavar="character"),
  
  make_option(c("-t", "--target_freq"), type="character", default="deniHCSS35unique", 
              help="name of the statistic that I want to Krig. E.g. deniHCSS35unique . ", metavar="character"),
  make_option(c("-r", "--resolution"), type="integer", default=500, 
              help="grid resolution at which to interpolate over the map, in km. The Kriging algorithm apparently scales as O(x^3) with x the number of grid points, so reducing this too far can really slow things down!", metavar="character"),
  make_option(c("-b", "--plot_buffer"), type="integer", default=500, 
              help="buffer around the most extreme datapoints (x/2 km in each direction)", metavar="character"),
  make_option(c("-j", "--jitter"), type="integer", default=10,
              help="jitter individual locations in a circle of this radius (km) to avoid identical locations. Thsi is default behaviour.", metavar="character"),
  make_option(c("-a", "--overlap_average"), type="logical", default=F, action="store_true",
              help="flag instructing the program to average individuals who are in the same location rather than jittering them.", metavar="character"),
  make_option(c("-g", "--log_krig"), type="logical", default=F, action="store_true",
              help="flag instructing the program to perform the kriging on the log10 target_freq. Can be used to detect subtle differences between regions.", metavar="character"),
  make_option(c("-c", "--continent"), type="character", default=NULL,
              help="instruction to plot the grid to a specific continent [europe, eurasia, asia, america, africa, isea, papua, world] rather than the adaptive map bounded around available datapoints.", metavar="character"),
  make_option(c("-p", "--populations_to_keep"), type="character", default=NULL,
              help="instruction to only include individuals in one of a list of populations. E.g. if you just want to draw New Guinea based on Papuan samples then you could add continent_papua or subset_papuaExclFrancois, or a list can be provided using e.g. continent_eisea,continent_papua.", metavar="character"),
  make_option(c("-v", "--validate_methods"), type="logical", default=F, action="store_true",
              help="flag instructing the program to subset and cross-validate data, printing the error in proximity ploygons, nearest neighbour, (fitted) inverse distance weighting and kriging methods.", metavar="character"),
  make_option(c("-d", "--datatype"), type="character", default="bp",
              help="data type - bp (divide by 10e6 to show MB; plot min to max); frequencies (don't change data; plot 0 to max); unchanged (don't change data; plot min to max).", metavar="character"),
  make_option(c("-m", "--nmax"), type="integer", default=-1,
              help="the number of nearest observations that should be used for a kriging prediction. A low number will radically speed up kriging and is very useful for testing, but may create anomalies.", metavar="character"),
  
  
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output file name for plotting.", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Select populations and ascertainment scheme. E.g. IDenDerFix_HumDerMax0_ADenAny_ANeaAnc.w5000.noAfPoly means (I think)  i) all introgressed haplotypes found in target population (New Guinea) for a genomic region to be fixed derived, ii) all human haplotypes found in target population (New Guinea) to be fixed ancestral, iii) the Altai Neanderthal to be fixed ancestral and iv) SNP isn't polymorphic in the African sample.
pops_to_compare <- c("deniHCSS35unique.subset_papuaExclFrancoisUVBaining.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc.w5000.noAfPoly", "deniHCSS35unique.subset_papuaExclFrancoisUVBaining.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc-resample0.w5000.noAfPoly", "deniHCSS35unique.continent_easia.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc.w5000.noAfPoly", "deniHCSS35unique.continent_easia.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc-resample0.w5000.noAfPoly", "deniHCSS35unique.continent_sasia.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc.w5000.noAfPoly", "deniHCSS35unique.continent_sasia.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc-resample0.w5000.noAfPoly", "neanHCSS35unique.subset_papuaExclFrancoisUVBaining.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc.w5000.noAfPoly", "neanHCSS35unique.subset_papuaExclFrancoisUVBaining.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc-resample0.w5000.noAfPoly", "neanHCSS35unique.continent_easia.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc.w5000.noAfPoly", "neanHCSS35unique.continent_easia.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc-resample0.w5000.noAfPoly", "neanHCSS35unique.continent_sasia.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc.w5000.noAfPoly", "neanHCSS35unique.continent_sasia.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc-resample0.w5000.noAfPoly", "neanHCSS35unique.continent_europe.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc.w5000.noAfPoly", "neanHCSS35unique.continent_europe.INeaDerFix_HumDerMax0_ANeaAny_ADenAnc-resample0.w5000.noAfPoly") # final ascertainment scheme
pops_to_compare_names <- c("Deni_Papua", "Deni_Papua_Resample", "Deni_EAsia", "Deni_EAsia_Resample", "Deni_SAsia", "Deni_SAsia_Resample", "Nean_Papua", "Nean_Papua_Resample", "Nean_EAsia", "Nean_EAsia_Resample", "Nean_SAsia", "Nean_SAsia_Resample", "Nean_Europe", "Nean_Europe_Resample")

# Initiate empty raster stack to store results
output_rasters <- stack()

# Withhold plotting
withhold_plots <- TRUE

for (i in 1:length(pops_to_compare)) {
  
  opt$location_file <- "/Users/hl636/Documents/Hirzi/Cambridge/Archaic sharing/Individual_GPS_locations_revised.txt"
  opt$target_freq <- pops_to_compare[i]
  #opt$target_freq <- "deniHCSS35unique.subset_papuaExclFrancoisUVBaining.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc.w5000.noAfPoly"
  #opt$target_freq <- "deniHCSS35unique.subset_papuaExclFrancoisUVBaining.IDenDerFix_HumDerMax0_ADenAny_ANeaAnc-resample0.w5000.noAfPoly" # random SNPs
  opt$freq_file <- "/Users/hl636/Documents/Hirzi/Cambridge/Archaic sharing/_combinedMeanFreqData_redoNoUVBaining.csv"
  opt$continent <- "oldworld"
  #opt$resolution <- 500
  opt$resolution <- 15
  opt$overlap_average <- F
  opt$frequencies <- T
  opt$datatype <- "frequencies"
  opt$jitter <- 10
  opt$plot_buffer <- 500
  # Note, populations_to_keep argument doesn't work as intended. Keeping this NULL for now and filtering for desired populations and individuals in a later step (just after reading in input files)
  #opt$populations_to_keep <- "continent_wisea,continent_eisea,subset_papuaExclFrancoisUVBaining,continent_oceania,continent_europe,continant_sasia,continent_easia,continent_africa"
  #opt$populations_to_keep <- NULL
  
  if (is.null(opt$continent) == F) {
    if (opt$continent %in% c("europe", "asia", "america", "africa", "eurasia", "isea", "papua", "world", "oldworld") == F) {
      stop("If continent is provided the implemented options are [europe, asia, america, africa, eurasia, isea, papua, world, oldworld]")
    }
  }
  
  if (opt$log_krig == T & opt$validate_methods == T) {
    stop("Cross-validation on logged data not implemented.")
  }
  
  ## Timer ##
  start_time <- Sys.time()
  
  # Function for working with read in locations and coverage data
  generate_loc_freq_by_ind = function( 
    locations_with_ind = list(), 
    freq_data_with_ind = list(),
    column = "",
    populations_to_keep_list = list()
  ){
    freq_data <- double()
    latitudes <- double()
    longitudes <- double()
    
    for (ind in locations_all$IND) {
      
      if (ind %in% freq_data_with_ind$Ind & (is.null(populations_to_keep_list) | length(intersect(populations_to_keep_list, unlist(strsplit(as.character(freq_data_all$Pops[freq_data_with_ind$Ind == ind]), "\\|")))) > 0)) {
        freq_obs <- freq_data_with_ind[[column]][freq_data_with_ind$Ind == ind]
        if (freq_obs != -1 & is.na(freq_obs) == F) {
          freq_data <- append(freq_data, freq_obs)
          latitudes <- append(latitudes, locations_all$LATITUDE[locations_all$IND == ind])
          longitudes <- append(longitudes, locations_all$LONGITUDE[locations_all$IND == ind])
        }
      }
      if (paste(ind,"_1",sep="") %in% freq_data_with_ind$Ind & (is.null(populations_to_keep_list) | length(intersect(populations_to_keep_list, unlist(strsplit(as.character(freq_data_all$Pops[freq_data_with_ind$Ind == paste(ind,"_1",sep="")]), "\\|")))) > 0)) {
        freq_obs <- freq_data_with_ind[[column]][freq_data_with_ind$Ind == paste(ind,"_1",sep="")]
        if (freq_obs != -1 & is.na(freq_obs) == F) {
          freq_data <- append(freq_data, freq_obs)
          latitudes <- append(latitudes, locations_all$LATITUDE[locations_all$IND == ind])
          longitudes <- append(longitudes, locations_all$LONGITUDE[locations_all$IND == ind])
        }
      }
      if (paste(ind,"_2",sep="") %in% freq_data_with_ind$Ind & (is.null(populations_to_keep_list) | length(intersect(populations_to_keep_list, unlist(strsplit(as.character(freq_data_all$Pops[freq_data_with_ind$Ind == paste(ind,"_2",sep="")]), "\\|")))) > 0)) {
        freq_obs <- freq_data_with_ind[[column]][freq_data_with_ind$Ind == paste(ind,"_2",sep="")]
        if (freq_obs != -1 & is.na(freq_obs) == F) {
          freq_data <- append(freq_data, freq_obs)
          latitudes <- append(latitudes, locations_all$LATITUDE[locations_all$IND == ind])
          longitudes <- append(longitudes, locations_all$LONGITUDE[locations_all$IND == ind])
        }
      }
    }
    locations = list("LATITUDE"=latitudes, "LONGITUDE" = longitudes)
    return(list("locations"=locations, "freq_data"=freq_data))
  }
  
  # Function for adding jitter in a cicle to lat/long locations in a dataframe. E.g. 10km in all directions.
  # This isn't very simple as i) I want a circle of e.g. 10kb radius and ii) lat/long degrees vary in km depending on how close to the equator you are
  add_jitter_lat_long = function( 
    locations = list(),
    radius = 10
  ){
    # https://gis.stackexchange.com/questions/25877/generating-random-locations-nearby
    latitudes_jittered = double()
    longitudes_jittered = double()
    for (ind_idx in 1:length(locations[,1])) {
      lat_long = locations[ind_idx,]
      r <- radius / 111.30 #Convert radius in km into degrees at the equator
      u <- runif(1,0,1)
      v <- runif(1,0,1)
      w <- r * sqrt(u)
      t <- 2 * pi * v
      x <- w * cos(t) 
      y <- w * sin(t)
      x_adj <- x / cos(lat_long[1] * (pi/180)) #Nb degrees to radians
      latitudes_jittered <- append(latitudes_jittered, as.double(lat_long[1] + y))
      longitudes_jittered <- append(longitudes_jittered, as.double(lat_long[2] + x_adj))
    }
    locations_jittered = list("LATITUDE"=latitudes_jittered, "LONGITUDE" = longitudes_jittered)
    return(as.data.frame(locations_jittered))
  }
  
  average_sameloc = function( 
    locations = list(),
    freq_data = list()
  ){
    # Go through each individual. If the lat/long hasn't been observed then add it in. If the lat/long has been observed, add in the freq and add 1 to the number of observations.
    lat_long <- locations[1,]
    latitudes <- c(as.double(lat_long[1]))
    longitudes <- c(as.double(lat_long[2]))
    freq_data_sum <- c(as.double(freq_data[1]))
    freq_data_num <- c(as.double(1))
    for (ind_idx in 2:length(locations[,1])) {
      lat_long = as.double(locations[ind_idx,])
      #print(lat_long)
      if (sum((lat_long[1] == latitudes) * (lat_long[2] == longitudes)) > 0) {
        # Observed before
        #print("Obs")
        obs_mask = (lat_long[1] == latitudes) * (lat_long[2] == longitudes) > 0
        freq_data_sum[obs_mask] = freq_data_sum[obs_mask] + as.double(freq_data[ind_idx])
        freq_data_num[obs_mask] = freq_data_num[obs_mask] + 1
      } else {
        # Not observed. Add to lat, long, freq_data_sum, freq_data_num
        #print("NotObs")
        latitudes = append(latitudes, as.double(lat_long[1]))
        longitudes = append(longitudes, as.double(lat_long[2]))
        freq_data_sum = append(freq_data_sum, as.double(freq_data[ind_idx]))
        freq_data_num = append(freq_data_num, 1)
      }
      #print(latitudes)
      #print(longitudes)
      #print(freq_data_sum)
      #print(freq_data_num)
    }
    #print(freq_data_sum)
    #print(freq_data_num)
    freq_data_av = freq_data_sum / freq_data_num
    locations_unique = list("LATITUDE"=latitudes, "LONGITUDE" = longitudes)
    return(list("locations"=as.data.frame(locations_unique), "freq_data"=freq_data_av))
    #return(NULL)
  }
  
  ### READ IN THE DATA ###
  
  # Applying to coverage data
  locations_all = read.table(opt$location_file, col.names = c("IND", "LATITUDE", "LONGITUDE"), check.names=FALSE)
  freq_data_all = read.table(opt$freq_file, sep = ",", header = T, check.names=FALSE)
  populations_to_keep_list = if(is.null(opt$populations_to_keep)) NULL else unlist(strsplit(opt$populations_to_keep, ","))
  
  # Since the populations_to_keep argument doesn't work as intended, we remove unwanted populations and individuals here
  #identical(freq_data_all$Ind, locations_all$IND)
  freq_data_all <- freq_data_all[freq_data_all$Pop != 'continent_america',]
  locations_all <- merge(freq_data_all, locations_all, by.x = "Ind", by.y = "IND")[,c(1, ncol(freq_data_all)+1, ncol(freq_data_all)+2)]
  colnames(locations_all)[1] <- "IND"
  #identical(freq_data_all$Ind, locations_all$IND)
  
  
  ### RETRIEVE THE REQUIRED DATA AND JITTER IF REQUESTED ###
  
  cat("Reading in location and frequency data for", opt$target_freq, "\n")
  
  loc_freq <- generate_loc_freq_by_ind(locations_all, freq_data_all, opt$target_freq, populations_to_keep_list)
  locations_unjittered <- as.data.frame(loc_freq$locations)
  if (opt$overlap_average == F) {
    locations <- add_jitter_lat_long(locations = locations_unjittered, radius = opt$jitter)
    freq_data <- loc_freq$freq_data
  } else {
    locations_freq_unique <- average_sameloc(locations = locations_unjittered, freq_data = loc_freq$freq_data)
    locations = locations_freq_unique$locations
    freq_data = locations_freq_unique$freq_data
  }
  
  if (opt$datatype == 'bp') {
    freq_data = freq_data / 1000000.0 #Conver to MB
  }
  
  cat(length(freq_data), "frequency/location pairs read in.\n")
  
  
  ### Prepare data formats required for different functions ###
  
  loc_freq_final = locations
  loc_freq_final$freq = freq_data
  loc_freq_final_sf <- st_as_sf(loc_freq_final, coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84 +units=km")
  loc_freq_final_sf = st_transform(loc_freq_final_sf, "+proj=merc +datum=WGS84 +units=km")
  
  
  ### Prepare grid for interpolation, and map for masking / plotting ###
  ext <- terra::ext(loc_freq_final_sf)
  
  ### BUILD THE WORLD GRID ###
  
  continent_list <- c(terra::ext(c(-20037.51, 20047.51, -8000, 14000)), #world
                      terra::ext(c(10000, 17000, -2250, 3000)), #isea
                      terra::ext(c(-1500, 5000, 3500, 12000)), #europe
                      terra::ext(c(3000, 20047, -1500, 12000)), #asia
                      terra::ext(c(-1500, 20047, 800, 12000)), #eurasia
                      terra::ext(c(-19000, -1700, -8000, 12000)), #america
                      terra::ext(c(-2000, 6700, -5000, 5000)), #africa
                      terra::ext(c(14500, 17000, -1500, 200)), #papua
                      terra::ext(c(-3000, 20047.51, -6000, 14000)) #old world
  )
  names(continent_list) <- c("world", "isea", "europe", "asia", "eurasia", "america", "africa", "papua", "oldworld")
  
  if (is.null(opt$continent) == T) {
    ext_buffer <- ext + opt$plot_buffer
  } else {
    ext_buffer <- continent_list[[opt$continent]]
  }
  
  # Use libarary rnaturalearth to build an Earth grid.
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sv")
  world_merc = terra::project(world, "+proj=merc +datum=WGS84 +units=km")
  
  x_range <- seq(ext_buffer[1], ext_buffer[2], by = opt$resolution)
  y_range <- seq(ext_buffer[3], ext_buffer[4], by = opt$resolution)
  grid <- expand.grid(x = x_range, y = y_range)
  grid_vect <- vect(grid, geom=c("x", "y"), crs = crs(loc_freq_final_sf))
  grid_sf <- st_as_sf(grid, coords = c("x", "y"))
  cat("There are", length(grid_sf$geometry) , "values in the kriging grid.\n")
  
  ### CALCULATE A VARIOGRAM ###
  
  # Use gstat to create an emperical variogram ‘v’
  
  pts <- vect(locations, geom=c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84 +units=km")
  pts_merc <- project(pts, "+proj=merc +datum=WGS84 +units=km")
  pts_merc_df <- data.frame(geom(pts_merc)[,c("x", "y")], as.data.frame(freq_data))
  
  gs <- gstat(formula=freq_data~1, locations = ~x+y, data = pts_merc_df)
  v <- variogram(gs, width=20)
  if (withhold_plots == FALSE) {
    plot(v)
  }
  
  
  ### DO THE AUTO KRIGING ###
  
  start_autoKrige = Sys.time()
  kp <- autoKrige(formula = if (opt$log_krig == F) {freq_data~1} else {log(freq_data)~1},
                  input_data = loc_freq_final_sf,
                  new_data = grid_sf,
                  #model = c("Exp", "Sph", "Gau", "Exc", "Mat", "Ste", "Cir", "Bes", "Pen", "Per", "Wav", "Hol", "Log", "Spl", "Leg"),
                  model = c("Exp", "Sph", "Gau", "Mat", "Ste"),
                  kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10),
                  start_vals = c(mean(v$gamma)/4, max(v$dist)/2, max(v$gamma)*0.9),
                  nmax = opt$nmax)
  
  end_autoKrige = Sys.time()
  cat("Kriging time in mins.", difftime(end_autoKrige, start_autoKrige, units = "mins"), "\n\n")
  
  kp$var_model
  if (withhold_plots == FALSE) {
    plot(variogramLine(kp$var_model, max(v[,2])*1.05), type='l', ylim=c(0,max(v[,3])*1.05))
    points(v[,2:3], pch=20, col='red')
  }
  
  
  ### CROSS-VALIDATION IF REQUESTED ###
  
  # This optionally cross-validates the kriging procedure against other options. Based on code in https://rspatial.org/raster/analysis/4-interpolation.html
  if (opt$validate_methods == T) {
    RMSE <- function(observed, predicted) {
      sqrt(mean((predicted - observed)^2, na.rm=TRUE))
    }
    null <- RMSE(mean(freq_data), freq_data) # Example prediction based on the average of data points
    
    # function for training the IDW fitting
    f1 <- function(x, test, train) {
      nmx <- x[1]
      idp <- x[2]
      if (nmx < 1) return(Inf)
      if (idp < .001) return(Inf)
      m <- gstat(formula=freq~1, locations=train, nmax=nmx, set=list(idp=idp))
      p <- predict(m, newdata=test, debug.level=0)$var1.pred
      RMSE(test$freq, p)
    }
    ndatapoints = length(loc_freq_final[[1]])
    ndatapoints_remainder5 = ndatapoints %% 5
    ndatapoints_quotient5 = ndatapoints %/% 5
    kf = c()
    for (k in 1:5) {
      if (ndatapoints_remainder5 >= k)
      {kf = append(kf, rep(k, ndatapoints_quotient5 + 1))
      }
      else {
        kf = append(kf, rep(k, ndatapoints_quotient5))
      }
    }
    rmsepp <- rep(NA, 5)
    rmsenn <- rep(NA, 5)
    rmseidw <- rep(NA, 5)
    rmsek <- rep(NA, 5)
    for (k in 1:5) {
      test <- loc_freq_final_sf[kf == k, ]
      train <- loc_freq_final_sf[kf != k, ]
      # proximity polygons
      v <- voronoi(vect(train))
      p <- extract(v, vect(test))
      rmsepp[k] <- RMSE(test$freq, p$freq)
      # neareset neighbour
      gscvnn <- gstat(formula=freq~1, locations=train, nmax=5, set=list(idp = 0))
      p <- predict(gscvnn, test)$var1.pred
      rmsenn[k] <- RMSE(test$freq, p)
      # inverse distance weighted, fitted
      opt_params <- optim(c(8, .5), f1, test=test, train=train)
      gscvidw <- gstat(formula=freq~1, locations=train, nmax=opt_params$par[1], set=list(idp=opt_params$par[2]))
      p <- predict(gscvidw, test)
      rmseidw[k] <- RMSE(test$freq, p$var1.pred)
      # kriging
      gscv <- gstat(formula=freq~1, locations=train)
      vcv <- variogram(gscv, width=20)
      kpcv <- autoKrige(formula = freq~1,
                        input_data = train,
                        new_data = test,
                        #model = c("Exp", "Sph", "Gau", "Exc", "Mat", "Ste", "Cir", "Bes", "Pen", "Per", "Wav", "Hol", "Log", "Spl", "Leg"),
                        model = c("Exp", "Sph", "Gau", "Mat", "Ste"),
                        kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10),
                        start_vals = c(mean(vcv$gamma)/4,max(vcv$dist)/2, max(vcv$gamma)*0.9))
      rmsek[k] <- RMSE(test$freq, kpcv$krige_output$var1.pred)
    }
    rmsepp_mi = 1 - mean(rmsepp)/null
    rmsenn_mi <- 1 - mean(rmsenn)/null
    rmseidw_mi <- 1 - mean(rmseidw)/null
    rmsek_mi <- 1 - mean(rmsek)/null
    print("Relative RMSE versus mean values for proximity polygons, 5-nearest neighbour, optimised IDW and kriging (higher is better):")
    c(rmsepp_mi, rmsenn_mi, rmseidw_mi, rmsek_mi)
    
    ### Plot the the optimised IDW fitting ###
    ndatapoints = length(loc_freq_final[[1]])
    ndatapoints_remainder2 = ndatapoints %% 2
    ndatapoints_quotient2 = ndatapoints %/% 2
    kf2 = c()
    for (k in 1:2) {
      if (ndatapoints_remainder2 >= k)
      {kf2 = append(kf2, rep(k, ndatapoints_quotient2 + 1))
      }
      else {
        kf2 = append(kf2, rep(k, ndatapoints_quotient2))
      }
    }
    opt_params <- optim(c(80, .5), f1, test=loc_freq_final_sf[kf2 == 1,], train=loc_freq_final_sf[kf2 == 2,])
    gsidw <- gstat(formula=freq_data~1, locations= ~x+y, data = pts_merc_df, nmax=opt_params$par[1], set=list(idp=opt_params$par[2]))
    idw <- interpolate(rast(grid_sf, nrows = length(y_range), ncols = length(x_range), crs = crs(loc_freq_final_sf)), gsidw)
    idw <- mask(idw, world_merc)
    if (withhold_plots == FALSE) {
      plot(idw)
    }
  }
  
  
  ### APPLYING KRIGING RESULTS TO THE MAP ###
  
  #These operations can actually take some time at fine resolution
  ok <- vect(kp$krige_output)
  ok <- mask(ok, world_merc)
  names(ok) <- c('prediction', 'variance', 'stdev')
  
  ok2 <- rast(ok, nrows = length(y_range), ncols = length(x_range))
  ok2 = rasterize(ok, y = ok2, field = "prediction", fun = "mean")
  names(ok2) <- 'prediction'
  
  
  ### SAVE THE OUTPUT ###
  if (is.null(opt$out) == F){
    # Open the saving file
    if (stri_sub(opt$out,-4,-1) == ".pdf") {
      pdf(opt$out)
    }
    else if (stri_sub(opt$out,-4,-1) == ".jpg") {
      jpeg(opt$out, width = 1000, height = 1000)
    }
    else if (stri_sub(opt$out,-4,-1) == ".png") {
      png(opt$out, width = 1080, height = 1000)
    }
    # Create the plot
    terra::plot(x = ok2, y = "prediction", col = map.pal("viridis", 100), type = "continuous", xlab = ("Kilometers from 0 longitude"), ylab = ("Kilometers from equator"), xlim = c(xmin(ext_buffer),xmax(ext_buffer)), ylim = c(ymin(ext_buffer), ymax(ext_buffer)), zlim = c(if (opt$datatype == "bp" | opt$datatype == "unchanged") {min(values(ok2), na.rm = T)} else 0.0, max(values(ok2), na.rm = T)))
    terra::plot(world_merc, add=TRUE)
    points(pts_merc[g], cex = .5, pch = 19, col = "red")
    # Close the output file
    invisible(dev.off())
  } else {
    # Just create the plot
    if (withhold_plots == FALSE) {
      terra::plot(x = ok2, y = "prediction", col = map.pal("viridis", 100), type = "continuous", xlab = ("Kilometers from 0 longitude"), ylab = ("Kilometers from equator"), xlim = c(xmin(ext_buffer),xmax(ext_buffer)), ylim = c(ymin(ext_buffer), ymax(ext_buffer)), zlim = c(if (opt$datatype == "bp" | opt$datatype == "unchanged") {min(values(ok2), na.rm = T)} else 0.0, max(values(ok2), na.rm = T)))
      terra::plot(world_merc, add=TRUE)
      points(pts_merc[g], cex = .5, pch = 19, col = "red")
    }
  }
  
  ## Save raster to stack
  ok3 <- raster(ok2)
  names(ok3) <- paste(pops_to_compare_names[i], names(ok3), sep="_")
  output_rasters <- raster::stack(output_rasters, ok3)
  
  ## Timer ##
  end_time <- Sys.time()
  
  cat("Done, total time in mins.", difftime(end_autoKrige, start_autoKrige, units = "mins"), "\n\n")
  
}

# Write out final results
#setwd("/Users/hl636/Documents/Hirzi/Cambridge/Archaic sharing/")
#writeRaster(output_rasters, "Archaic_maps_raw_revised")

