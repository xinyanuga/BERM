###This script applies already created aboveground biophysical models to new 
###prediction data sets, usually broad scale site level data
library(R.utils); library(tidyverse);library(RColorBrewer)
library(data.table); library(xgboost); library(randomForest)
library("paradox");  library("mlr3");library("mlr3learners"); 
library("mlr3tuning"); library("randomizr")


#####Site level prediction set data####
###Read in processed Landsat 8 data for the site level prediction set ####
newpixels<-read_csv("output/gce_pixels_landsat8_processed.csv")
newpixels$date<-as.Date(newpixels$date)

##Read in the greenup and surface temperature estimates for the site level prediction set
newpixels.heat.l8<-read_csv("output/gce_pixels_landsat8_greenup_via_lst_and_air_lags.csv")
newpixels.heat.l8$greenup<-as.Date(newpixels.heat.l8$greenup.soil)

##variable names for use in model prediction for aboveground biophysical models
feature.names<-c("b1","b2","b3","b4","b5","b6","b7","ndvi", "pheno", "vari")

##note for each of the models below, the following column names are expected as predictor features:
#feature.names<-c("b1","b2","b3","b4","b5","b6","b7","ndvi", "pheno", "vari")
##the above was handled by the newpixels processing script

####Predict aboveground biomass####
##first load the AGB model
agmod1<-xgb.load("output/xgb_agb_model")

##now predict the data
newpixels<-data.frame(newpixels)
newpixels$predag.allom.l8 <-(predict(agmod1, data.matrix(newpixels[,feature.names])))

####Predict foliar N####
##first load the foliar N model
Nmod1<-xgb.load("output/xgb_foliarN_model")
newpixels<-data.frame(newpixels)

##now predict the data
newpixels$predn.l8<- (predict(Nmod1, data.matrix(newpixels[,feature.names])))


####Predict LAI####
##first load the model
laimod1<-xgb.load("output/xgb_lai_model")

##now predict the data
newpixels<-data.frame(newpixels)
newpixels$predlai.l8<- (predict(laimod1, data.matrix(newpixels[,feature.names])))


####Predict CHL####
##first load the model
feature.names<-c("b1","b2", "b3","b4","b5","b6","b7", "gari", "ndmi", "pheno", "vari", "ndvi")
chlmod1<-xgb.load("output/xgb_chl_model")

##now predict the data
newpixels<-data.frame(newpixels)
newpixels$predchl.l8<- (predict(chlmod1, data.matrix(newpixels[,feature.names])))


##save a copy of the prediction results for later use
write_csv(newpixels, path = "output/xgb_predicted_biophysical_seagrant_landsat8_allpixels.csv")
