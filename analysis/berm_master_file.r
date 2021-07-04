####Belowground Ecosystem Resilency Model (BERM) Master file###
###Run code in the order below to generate Belowground biomass estimates for Landsat 8 pixels known to represent Spartina alterniflora marsh

###change to your project working directory
### This assumes the following folders will be present:
###  analysis (for R code), data (for raw data), functions (for custom function scripts)
###  output (for data processed via R scripts or generated via R models),  results (for plots/visuals)

setwd("/home/jessica/UGA/Projects/BERM/")

##load in needed custom functions and models
source("functions/soil_temp_from_lst_and_air_lags.r")
source("../functions/calc_lst_l8.r")
source("../functions/calc_indices_l8.r")
source ('../functions/colorRampPaletteAlpha.r')
source("analysis/spad_to_chl.r")

##tidal flooding prediction model
super_model <- readRDS("output/random_forest_flood_model_l8.rds")

##This script retrieves atmospheric correction parameters for thermal band data in order to calculate land surface temperature
##note that this script contacts the NASA webtool and retrieves data, so don't run this frivolously, uncomment the line below to run
##data values apply to a broad spatial scale (1 degree of lat/long), so we don't have to get this for every pixel, only if we switch locations to a new lat/long
#source("analysis/get_atmospheric_correction.r")


##Take the raw vegetation field data and interpolate them so that we can merge with the nearest cloud and flood free Landsat observation date
source("analysis/smooth_veg_timeseries.r")

##Convert daily gridded Daymet data, acquired from GEE, to monthly averages
source("analysis/process_daymet.r")

###create and save xgboost models for aboveground proxies 
source("analysis/xgboost_landsat_biophysical.r")

###to predict aboveground biophysical models on new data, usually site level data
##first process the new data so the right variables are present
source("analysis/process_site_level_newpixels.r")

##Next predict the biophysical data on the new pixels
source("analysis/predict_biophysical_data_newpixels.r")

##Next prepare the output of the xgboost landsat biophyical script to build the 
##belowground biomass model
source("analysis/process_groundtruth_dataset_for_xgboost_bgbiomass_model.r")

##Next prepare the site level newpixel data that don't have ground-truth info  
## to predict from belowground biomass model
source("analysis/process_newpixels_for_xgb.r")

##Next create the belowground biomass model through a spatially crossvalidated 
##nested resampling scheme
##note code for predicting the site level novel data is at the bottom of this script
source("analysis/model_xgb_belowground_biomass.r")

##If desired, convert the output of the above script into raster files and then 
##time series maps
source("analysis/create_raster_timeseries.r")

