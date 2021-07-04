###This script creates the aboveground biophysical models from the ground-truth field data
###It also estimates pixel flood status and percent flooded, date of pixel greenup and LST calculations

library(R.utils); library(caret); library(tidyverse);library(RColorBrewer)
library(data.table); library(xgboost); library(randomForest)
library("paradox");  library("mlr3");library("mlr3learners"); 
library("mlr3tuning"); library("randomizr")

##output and save new models and data sets?? Set to True/False below. 
##Note that True will overwrite the existing models and output datasets, so you'll want a backup if you don't intend this
save_model<-FALSE
#save_model<-TRUE

color<-brewer.pal(8, "Dark2")
trans.color<-addalpha(color,0.9)


####Ground-truth data for field observations####
##Read in location data for vegetation ground truth plots
plots<-read.csv("data/biomass_plots.csv", header=T)
##not all dates have all pixels, so we need to make the averages here in order to avoid duplicate elevations later
pixels<- plots %>% group_by(pix) %>% summarise(elevation=mean(elevation), 
                                               utm_east=mean(utm_east), 
                                               utm_north=mean(utm_north), 
                                               lat=mean(lat), long=mean(long))
plots<-dplyr::select(plots, plot, pix)

###Read in Landsat8 observations for ground-truth plots, collected from Google Earth Engine Script####
### read in the landsat data for the main vegetation data plots
landsat<-read.csv("data/L8_seagrant_plots_2013_2020_noclouds.csv",
                  header = T)
### read in Landsat data for areas where we had the temperature data plots
lsat<-read.csv("data/L8_tidbit_plots_20170108_noclouds.csv", stringsAsFactors = F)
lsat$long<-substr(lsat$.geo, 32,49)
lsat$lat<-substr(lsat$.geo, 51,67)
lsat<-dplyr::select(lsat, names(lsat)[names(lsat) %in% names(landsat)])
lsat$plot<-paste0("tb", lsat$plot)
##combine the temperature and vegetation data into a single data frame
landsat<-plyr::rbind.fill(landsat, lsat)

##process the landsat data to create needed time, date, and location columns
landsat$time<-substr(landsat$date,12,19)
landsat$date<-as.Date(substr(landsat$date,1,10))
landsat$year<-as.numeric(format(landsat$date, "%Y"))
landsat$doy<-as.numeric(format(landsat$date, "%j"))
landsat$mo<-as.numeric(format(landsat$date, "%m"))
landsat<-merge(landsat, plots, by="plot")
##filter the landsat 8 data if needed through the pixel_qa landsat mask; 
##note that the current Landsat 8 Google Earth Engine script only returns cloud filtered data
##see https://landsat.usgs.gov/landsat-surface-reflectance-quality-assessment
landsat$pixel_qa<-intToBin(landsat$pixel_qa)
##"000010" from the right means no fill, yes clear, no water, no cloud shadow
## no water above is via Landsat 8's pixel qa mask, and misses a lot of marsh flooding in mixed pixels, we'll handle this below
landsat$qa_good<-ifelse(str_sub(landsat$pixel_qa,-6,-1)=="000010",T,F) 

## subset to only good pixels as indicated by the pixel_qa mask
landsat<-landsat[landsat$qa_good==T,]

##remove columns we don't need now from Landsat 8 Google Earth Engine data
landsat<-dplyr::select(landsat,-c(system.index, .geo, radsat_qa, sr_aerosol, pixel_qa, qa_good))
landsat$year<-as.numeric(format(landsat$date, "%Y"))

## Scale the landsat 8 bands to represent surface reflectance
landsat$b1<-landsat$B1*0.0001
landsat$b2<-landsat$B2*0.0001
landsat$b3<-landsat$B3*0.0001
landsat$b4<-landsat$B4*0.0001
landsat$b5<-landsat$B5*0.0001
landsat$b6<-landsat$B6*0.0001
landsat$b7<-landsat$B7*0.0001

##this is the needed scalar for these thermal bands: B10 and B11
landsat$b10<-landsat$B10*0.1
landsat$b11<-landsat$B11*0.1

##clean up columns to remove the ones we now don't need
landsat<-dplyr::select(landsat, -c(B1,B2,B3,B4,B5,B6,B7,B10,B11))

##negative reflectance values can occur at scene edges and should be removed
landsat$b1<-ifelse(landsat$b1<0,NA, landsat$b1)
landsat$b2<-ifelse(landsat$b2<0,NA, landsat$b2)
landsat$b3<-ifelse(landsat$b3<0,NA, landsat$b3)
landsat$b4<-ifelse(landsat$b4<0,NA, landsat$b4)
landsat$b5<-ifelse(landsat$b5<0,NA, landsat$b5)
landsat$b6<-ifelse(landsat$b6<0,NA, landsat$b6)
landsat$b7<-ifelse(landsat$b7<0,NA, landsat$b7)

##filter to non-NA values
landsat<-landsat[is.na(landsat$b1)==F,]; landsat<-landsat[is.na(landsat$b7)==F,]; 
landsat<-landsat[is.na(landsat$b5)==F,];landsat<-landsat[is.na(landsat$b6)==F,]

##add in a site variable for grouping the data
landsat$site<-"fluxa"
landsat$site<-ifelse(str_sub(landsat$plot,1,2)=="ug", "ugami", landsat$site)
landsat$site<-ifelse(str_sub(landsat$plot,1,2)=="sk", "skida", landsat$site)
landsat$site<-ifelse(str_sub(landsat$plot,1,2)=="fb", "fluxb", landsat$site)
landsat$site<-ifelse(str_sub(landsat$plot,1,2)=="tb", "xtidb", landsat$site)

plot(landsat$date, landsat$b3)

###calculate a set of standard vegetation and spectral reflectance indices for landsat 8 data using our custom function
indices<-calc_index_l8(landsat)
landsat<-cbind(landsat, indices)

##read in variables needed to atmospherically correct thermal data; 
##These were generated by the get_atmospheric_correction.r script
atm<-read.csv("output/atmospheric_correction_varsL8_gce.csv")
atm$date<-as.Date(atm$date)
landsat<-left_join(landsat, atm, by="date")

##calculate lst with custom function
landsat$lst<-calc_lst(landsat)
landsat<-landsat[!is.na(landsat$lst),]

##create rolling pheno mean
##trick to get all dates, just create the full dates dataframe and merge it to the missing dates dataframe
y<-data.frame(date=seq.Date(from=min(landsat$date, na.rm=T), to=max(landsat$date, na.rm=T), by=1))
y<-expand.grid(date=y$date, plot=unique(landsat$plot))
landsat<-full_join(landsat,y, all.y=TRUE)
landsat<-arrange(landsat, plot,date)
##this is an even window roll, because there's an observation for every day
landsat<-landsat %>% dplyr::group_by(plot) %>%
  dplyr::mutate(pheno2= RcppRoll::roll_mean(pheno, n=30, na.rm=T,align="center", fill=NA))
landsat<-landsat[!is.na(landsat$b2),]

##Use the tidal flooding prediction random forest model (super_model); 
##We'll predict tidal flooding and filter the data to dry observations
print(super_model)

## create some more of the needed columns the flood model expects
landsat$mo<-factor(landsat$mo, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))

##predict flooding and then process the prediction column to interpretable labels
landsat$flooded <- predict(super_model, newdata=landsat)
landsat$est.flood<-ifelse(landsat$flooded==0, "dry", "wet")
landsat$est.flood<-factor(landsat$est.flood)
table(landsat$date, landsat$flooded)

##read in the smoothed vegetation data, to merge with the landsat data
data<-read.csv("output/smoothed_veg_timeseries_weekly.csv", header=T, stringsAsFactors = F)
data<-arrange(data, date)
data$date<-as.Date(data$date)
data$field_lai<-data$LAI
data<-dplyr::select(data, -year, -LAI)
data$site<-ifelse(str_sub(data$plot,1,2)=="tb", "xtidb", data$site)

##quality control on pixels
x<-table(landsat$date, landsat$est.flood);
##this gets rid of dates with many wet obs
x<-x[x[,2]<10,] 
##this last gets rid of dates with few dry satelitte obs.
x<-x[x[,1]>5,] 
landsat<-landsat[landsat$date %in% as.Date(row.names(x)),]
landsat<-landsat[landsat$est.flood=="dry",]

landsat<-landsat[is.na(landsat$b5)==F,]
landsat<-arrange(landsat, date)

#create a better heat interpolation
heat.l8<-setup_heat(landsat[landsat$site!="xtidb",], gce=F)
if(save_model==T){
  write_csv(heat.l8, "output/seagrant_pixels_landsat8_heat_processed.csv")
}


###Merge lsat and field landsat < 10 days apart
##find the week in landsat landsat that is closest to a week in field survey  to create a merge key
data$start<-data$date-10; data$end<-data$date+11
data<-data.table(data)
setkey(data, plot, start, end)

##add on empty rows with the date and plot of early landsat observations so that these make it into the merged landsat
earlydates<-unique(landsat[,c("date", "plot")])
earlydates<-earlydates[earlydates$date<min(data$start),]
var<-names(data)
tmp<-data[1,]
tmp2<-tmp
for (i in 2:nrow(earlydates)){
  tmp<-rbind(tmp, tmp2)
}
tmp[,1:ncol(tmp)]<-NA
tmp$date<-earlydates$date
tmp$plot<-earlydates$plot
data<-rbind(data, tmp)

landsat<-landsat[!(landsat$plot %in% c("tb1", "tb2")),]
landsat<-dplyr::select(landsat, -mo)
data<-dplyr::select(data, -mo, -site)
data<-data[!is.na(data$date),]

data$start<-as.Date(data$date)-10; data$end<-as.Date(data$date)+11
data<-data.table(data)
setkey(data, plot, start, end)

##now create an index in the lsat landsat, here we'll just repeat the date column
##twice so that if the landsat date is within the field interval, we'll get a match
landsat$lsat_date<-as.Date(landsat$date)
landsat$start<-landsat$date; landsat$end<-landsat$date
landsat<-dplyr::select(landsat, -date)

###this removes dates before the field landsat; but we want to predict those dates to get rolling means
#landsat<-landsat[landsat$lsat_date>=min(data$start),]
landsat<-data.table(landsat)
setkey(landsat, plot, start, end)


#join if the date of landsat is within the end and start interval of field,
##return nothing for rows without matches, return two rows for multiple matches (we'll pick the best one next)
data<-foverlaps(landsat, data, type="any", by.x=c("plot", "start", "end"),
                by.y=c("plot", "start", "end"), nomatch=NA, mult="all")

#return veg data and landsat data back to tibbles;
#data.table objects have some different behavior than data.frames and tibbles and require care when using
landsat<-as_tibble(landsat); data<-as_tibble(data)

##subset field obs to those that are close to a landsat date,
##elemintates duplicate field estimates for the same landsat date
data$date<-as.Date(ifelse(is.na(data$date),data$lsat_date, data$date), origin="1970-01-01")
data$diff<-abs(data$date-data$lsat_date)

data<-data %>% dplyr::group_by(lsat_date, plot) %>%
  filter(diff==min(diff, na.rm=T))
##we don't need this, already decided above the allowed date range
#data<-data %>%  filter(abs(diff)< 10)
data<-data %>% dplyr::select(-c(start,end, i.start, i.end, diff))

data$landsat_doy<-data$doy
data$mo<-as.numeric(format(data$date, "%m"))

##create pixel averages and pixel average flood status
data<-dplyr::mutate(data, pixel=paste(year,landsat_doy,pix, site,sep="_"))
data<-ungroup(data)

##summarize flood status, to write out freq of flood observations by pixel for later use
flood<-data %>%
  dplyr::group_by(date, pix,site, flooded) %>% dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)
flood$flooded<-as.numeric(flood$flooded)
flood$flooded<-ifelse(flood$flooded==1,0,flood$flooded)
flood$flooded<-ifelse(flood$flooded==2,1,flood$flooded)
flood<- flood %>% dplyr::group_by(pix) %>% dplyr::summarise(flood_time= sum(as.numeric(flooded), na.rm=T)
                                                            /length(flooded[!is.na(flooded)]))
if(save_model==T){
  write_csv(flood, "output/flooded_percent_from_landsat_seagrant_pixels.csv")
}

##get rid of variables we don't need now to prepare for pixel averaging of only dry pixels
data<-data %>% dplyr::select(-c( plot,est.flood,  time, zone, flooded,pixel,mo, 
                                 year,lsat_date,landsat_doy, elevation, utm_east, utm_north, lat, long
                                                                  )) 
data<-data %>% dplyr::group_by(date, pix,site) %>% dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)
data<-ungroup(data)
data$mo<-as.numeric(format(data$date, "%m"))
data$year<-as.numeric(format(data$date, "%Y"))
data$doy<-as.numeric(format(data$date, "%j"))
data$group<-paste(data$site, data$mo)
data<-data[!is.na(data$site),]
data<-data[!is.na(data$b1),]
data$id<-as.numeric(row.names(data))
data$spad<-data$chlorophyll
data$chlorophyll<-(predict(spad, newdata=data))

##XGBoost
##Set up the resampling scheme for nested resampling
outer<-1
inner<-3

clusters <-data$group
blocks <-  data$site
rep<-1:outer
vars<-paste0("cv",rep)
cvs<-matrix(data=NA,nrow=nrow(data),ncol=length(vars))

###create the cross valdidated folds####
###999 
set.seed(29382) 
for( i in 1: outer){
  cvs[,i]<- block_and_cluster_ra(blocks = blocks,
                                 clusters = clusters,
                                 prob_each = c(.35, .65))
}       

##require tidbit lai to go into the training data, as we only have one month for that location
cvs<-data.frame(cvs); names(cvs)<-vars
head(table(clusters, cvs$cv1));head(table(blocks, cvs$cv1))

#clusters; blocks
table(cvs$cv1)
head(table(blocks, cvs$cv1))

##set up training/testing index for landsat data for later access
for (i in 1:outer){
  data[[paste0("train.", i)]] <- cvs[,i]}

table(data$train.1)
table(cvs$cv1)


##create a list of parameters to focus on tuning, must specify type (Int=Integer, Dbl=continuous number, Fct=factor/character)
tune_ps = ParamSet$new(list(
  ParamDbl$new("eta", lower = 0.05, upper = 0.2),
  ParamInt$new("max_depth", lower = 4, upper = 11),
  ParamInt$new("min_child_weight", lower = 1, upper = 10),
  ParamDbl$new("subsample", lower = 0.25, upper = 0.8),
  ParamDbl$new("colsample_bytree", lower = 0.25, upper = 0.8)
))
tune_ps

###AGBIOMASS#####
tr<-data
tr<-dplyr::filter(tr, agm2.allom>10)
feature.names<-c("b1","b2","b3","b4","b5","b6","b7","ndvi", "pheno", "vari")
resp<-"agm2.allom"
tr<-data[data$train.1==1,c(resp, feature.names)]
tr<-tr[complete.cases(tr),]

traintask = TaskRegr$new(id = "agbiomass", backend = tr[,c(resp, feature.names)], target = resp)
traintask$select(feature.names)

##set up learner
##select a resampling stategy for the datasets and preformance measure to use to tune the paramters
##set up how the tuner should iterate through the hyperparamters to try; 
##set up how the tuning should stop: here it stops after 20 iterations
##see other options at: https://mlr3book.mlr-org.com/tuning.html
learner = lrn("regr.xgboost", predict_type = "response",  booster="gbtree",  objective   = "reg:squarederror", nrounds=500)
#resamp=rsmp("holdout", ratio=0.85)
resamp=rsmp("cv", folds=inner)
measure = msr("regr.mse")
tuner = mlr3tuning::tnr("grid_search", resolution = 25)
terminator = trm("stagnation", iters=4, threshold=0.1)

instance = mlr3tuning::TuningInstanceSingleCrit$new(
  task = traintask,
  learner = learner,
  resampling = resamp,
  measure = measure,
  search_space = tune_ps,
  terminator = terminator
)
instance

##run the optimization
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
#instance$result_learner_param_vals

##store the best parameters back in the learner
learner$param_set$values = instance$result_learner_param_vals
learner$train(traintask)

##see the result, eg, the best configuration of parameters found during the iterations
param <- instance$result_x_domain

resp<-tr$agm2.allom
clf <- xgboost(data        = data.matrix(tr[,c(feature.names)]),
               label       = resp,
               booster="gbtree",
               nrounds     = 400,
               params = param,
               verbose=0)

mat <- xgb.importance (feature_names = feature.names,model = clf)
xgb.plot.importance (importance_matrix = mat[1:18]) 
p <-predict(clf, data.matrix(tr[, feature.names]))
plot(resp, p, ylim=c(0,max(c(p, resp), na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T))); abline(0,1)
sqrt((sum((resp-p)^2, na.rm=T))/length(p))

resp<-data$agm2.allom[data$train.1==0]
p <-predict(clf, data.matrix(data[data$train.1==0, feature.names]))
rmse<-sqrt((sum((resp-p)^2, na.rm=T))/length(p))
print(rmse)
plot(resp, p, ylim=c(0,max(c(p, resp),na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T)), pch=19,
     col=trans.color[c(8,1,2,4)][factor(data$site[data$train.1==0&data$site!="xtidb"])],
     cex.lab=1.5, cex.axis=1.2)
abline(0,1)
mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.1f",round(rmse,1))) ~ g ~ m^-2), cex=1.2)
legend("bottomright", legend=unique(data$site[!is.na(data$site)&data$site!="xtidb"]), 
       title=expression(underline(site)),title.adj = 0.25,
       col=trans.color[c(8,1,2,4)][factor(unique(data$site[!is.na(data$site)&data$site!="xtidb"]))], pch=15, bty="n")

##aboveground biophysical results
#rmse[n]
sqrt((sum((resp-p)^2, na.rm=T))/length(resp))
#nRMSE
(sqrt((sum((resp-p)^2, na.rm=T))/length(resp)))/(max(resp, na.rm=T)-min(resp, na.rm=T))
#COV RMSE (alternate nRMSE)
(sqrt((sum((resp-p)^2, na.rm=T))/length(resp)))/(mean(resp, na.rm=T))


#mae[n]
(sum((p-resp), na.rm=T))/length(resp[!is.na(resp)])
#cor
cor(resp, p, use="pairwise.complete.obs")
#N
nrow(data[!is.na(data$agm2.allom)&data$train.1==0,])
nrow(data[!is.na(data$agm2.allom)&data$train.1==1,])

agmod1<-clf 
#clf<-agmod1
data<-data.frame(data)
data$predag.allom.l8<-(predict(agmod1, data.matrix(data[,feature.names])))
maxx<-max(c(data$agm2.allom, data$predag.allom.l8), na.rm=T)
minn<-min(c(data$agm2.allom, data$predag.allom.l8), na.rm=T)
sqrt((sum((data$agm2.allom-data$predag.allom.l8)^2, na.rm=T))/length(data$predag.allom.l8[!is.na(data$predag.allom.l8)]))

if(save_model==TRUE){
  ##save a copy of the model for future use
  saveRDS(agmod1, "output/xgb_agb_model.rds")
  xgb.save(agmod1,  "output/xgb_agb_model")
  agmod1<-xgb.load("output/xgb_agb_model")
  
  ##create a graph of measured vs predicted for training and testing data
  jpeg(file="results/xgb_l8_ag_allometeric.jpg",
       height=5, width=5, res=200, units="in")
  par(mar=c(4,5,0.5,0.75), oma=c(0,0.2,0.5,0.5), mfrow=c(1,1))
  plot(data$agm2.allom, data$predag.allom.l8,
       xlab=expression(paste("measured AG biomass g ", m^-2, sep="")),
       ylab =expression(paste("predicted AG biomass g ", m^-2, sep="")) ,
       ylim=c(minn,maxx),xlim=c(minn,maxx),  pch=c(1,19)[factor(data$train)],
       cex=c(1,0.8)[factor(data$train)], col=trans.color[c(8,1,2,4)][factor(data$site)],
       cex.lab=1.5, cex.axis=1.2)
  abline(0,1)
  mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(round(rmse,1)) ~ g~m^-2), cex=1.2)
  legend("bottomright",inset=c(0,0.07),legend=c("train", "test"), pch=c(19,1), bty="n", title=expression(underline(data)), title.adj = 0.3)
  legend("bottomright", inset=c(0.2,0),legend=unique(data$site[!is.na(data$site)&data$site!="xtidb"]),
         title=expression(underline(site)),title.adj = 0.25,
         col=trans.color[c(8,1,2,4)][factor(unique(data$site[!is.na(data$site)&data$site!="xtidb"]))], pch=15, bty="n")
  dev.off()
}



###PERCENT N#####
###for some reason this doesn't create a prediction with the new MLR3 package, so fitting model with old MLR package
tr<-data
tr<-dplyr::filter(tr, agm2.allom>10)
feature.names<-c("b1","b2","b3","b4","b5","b6","b7","ndvi", "pheno", "vari")
resp<-"percentN"
tr<-data[data$train.1==1,c(resp, feature.names)]
tr<-tr[complete.cases(tr),]

traintask = TaskRegr$new(id = "id", backend = tr[,c(resp, feature.names)], target = resp)
traintask$select(feature.names)

##set up learner
learner = lrn("regr.xgboost", predict_type = "response",  booster="gbtree",  objective   = "reg:squarederror", nrounds=500)
#resamp=rsmp("holdout", ratio=0.85)
resamp=rsmp("cv", folds=inner)
measure = msr("regr.mse")
tuner = mlr3tuning::tnr("grid_search", resolution = 25)
terminator = trm("stagnation", iters=4, threshold=0.1)

instance = mlr3tuning::TuningInstanceSingleCrit$new(
  task = traintask,
  learner = learner,
  resampling = resamp,
  measure = measure,
  search_space = tune_ps,
  terminator = terminator
)
instance

##run the optimization
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
#instance$result_learner_param_vals

##store the best parameters back in the learner
learner$param_set$values = instance$result_learner_param_vals
learner$train(traintask)

##see the result, eg, the best configuration of parameters found during the iterations
param <- instance$result_x_domain

resp<-tr$percentN
clf <- xgboost(data        = data.matrix(tr[,c(feature.names)]),
               label       = resp,
               booster="gbtree",
               nrounds     = 400,
               params = param,
               verbose=0)

mat <- xgb.importance (feature_names = feature.names,model = clf)
xgb.plot.importance (importance_matrix = mat[1:18]) 
p <-predict(clf, data.matrix(tr[, feature.names]))
plot(resp, p, ylim=c(0,max(c(p, resp), na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T))); abline(0,1)
sqrt((sum((resp-p)^2, na.rm=T))/length(p))

resp<-data$percentN[data$train.1==0]
p <-predict(clf, data.matrix(data[data$train.1==0, feature.names]))
rmse<-sqrt((sum((resp-p)^2, na.rm=T))/length(p))
print(rmse)
plot(resp, p, ylim=c(0,max(c(p, resp),na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T)), pch=19,
     col=trans.color[c(8,1,2,4)][factor(data$site[data$train.1==0&data$site!="xtidb"])],
     cex.lab=1.5, cex.axis=1.2)
abline(0,1)
mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.1f",round(rmse,1))) ~ g ~ m^-2), cex=1.2)
#legend("bottomright",inset=c(0,0.1),legend=c("train", "test"), pch=c(19,1), bty="n", title=expression(underline(data)), title.adj = 0.3)
legend("bottomright", legend=unique(data$site[!is.na(data$site)&data$site!="xtidb"]), 
       title=expression(underline(site)),title.adj = 0.25,
       col=trans.color[c(8,1,2,4)][factor(unique(data$site[!is.na(data$site)&data$site!="xtidb"]))], pch=15, bty="n")

##aboveground biophysical results
#mae[n]
(sum((p-resp), na.rm=T))/length(resp[!is.na(resp)])
##RMSE
sqrt((sum((resp-p)^2, na.rm=T))/length(resp))
#nRMSE
(sqrt((sum((resp-p)^2, na.rm=T))/length(resp)))/(max(resp, na.rm=T)-min(resp, na.rm=T))
#COV RMSE (alternate nRMSE)
(sqrt((sum((resp-p)^2, na.rm=T))/length(resp)))/(mean(resp, na.rm=T))
#cor
cor(resp, p, use="pairwise.complete.obs")
#Ntest
nrow(data[!is.na(data$percentN)&data$train.1==0,])
#Ntrain
nrow(data[!is.na(data$percentN)&data$train.1==1,])


Nmod1<-clf 
#clf<-Nmod1
data<-data.frame(data)
data$predn.l8<-(predict(Nmod1, data.matrix(data[,feature.names])))
maxx<-max(c(data$percentN, data$predn.l8), na.rm=T)
minn<-min(c(data$percentN, data$predn.l8), na.rm=T)
sqrt((sum((data$percentN-data$predn.l8)^2, na.rm=T))/length(data$predn.l8[!is.na(data$predn.l8)]))


if(save_model==TRUE){
  ##save a copy of the model for future use
  saveRDS(Nmod1, "output/xgb_foliarN_model.rds")
  xgb.save(Nmod1,  "output/xgb_foliarN_model")
  Nmod1<-xgb.load("output/xgb_foliarN_model")
  
  ##create a graph of measured vs predicted for training and testing data
  jpeg(file="results/xgb_l8_N.jpg", res=200, units="in",
     height=5, width=5)
  par(mar=c(4,5,0.5,0.75), oma=c(0,0.2,0.5,0.5), mfrow=c(1,1))
  plot(data$percentN, data$predn.l8, 
     xlab="measured % foliar N", 
     ylab ="predicted % foliar N",
     ylim=c(0,maxx),xlim=c(0,maxx),  pch=c(1,19)[factor(data$train)], 
     cex=c(1,0.8)[factor(data$train)], col=trans.color[c(8,1,2,4)][factor(data$site)],
     cex.lab=1.5, cex.axis=1.2)
  abline(0,1)
  mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.1f",round(rmse,1))) ~ "%"), cex=1.2)
  legend("bottomright",inset=c(0,0.07),legend=c("train", "test"), pch=c(19,1), bty="n", title=expression(underline(data)), title.adj = 0.3)
  legend("bottomright", inset=c(0.2,0),legend=unique(data$site[!is.na(data$site)&data$site!="xtidb"]), 
       title=expression(underline(site)),title.adj = 0.25,
       col=trans.color[c(8,1,2,4)][factor(unique(data$site[!is.na(data$site)&data$site!="xtidb"]))], pch=15, bty="n")
  dev.off()
}


####LAI#####
tr<-data
tr<-dplyr::filter(tr, agm2.allom>10)
feature.names<-c("b1","b2","b3","b4","b5","b6","b7", "ndvi", "pheno", "vari")
resp<-"field_lai"
tr<-data[data$train.1==1,c(resp, feature.names)]
tr<-tr[complete.cases(tr),]

traintask = TaskRegr$new(id = "id", backend = tr[,c(resp, feature.names)], target = resp)
traintask$select(feature.names)

##set up learner
learner = lrn("regr.xgboost", predict_type = "response",  booster="gbtree",  objective   = "reg:squarederror", nrounds=500)
#resamp=rsmp("holdout", ratio=0.85)
resamp=rsmp("cv", folds=inner)
measure = msr("regr.mse")
tuner = mlr3tuning::tnr("grid_search", resolution = 25)
terminator = trm("stagnation", iters=4, threshold=0.1)

instance = mlr3tuning::TuningInstanceSingleCrit$new(
  task = traintask,
  learner = learner,
  resampling = resamp,
  measure = measure,
  search_space = tune_ps,
  terminator = terminator
)
instance

##run the optimization
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
#instance$result_learner_param_vals

##store the best parameters back in the learner
learner$param_set$values = instance$result_learner_param_vals
learner$train(traintask)

##see the result, eg, the best configuration of parameters found during the iterations
param <- instance$result_x_domain

resp<-tr$field_lai
clf <- xgboost(data        = data.matrix(tr[,c(feature.names)]),
               label       = resp,
               booster="gbtree",
               nrounds     = 100,
               params = param,
               verbose=0)

mat <- xgb.importance (feature_names = feature.names,model = clf)
xgb.plot.importance (importance_matrix = mat[1:18]) 
p <-predict(clf, data.matrix(tr[, feature.names]))
plot(resp, p, ylim=c(0,max(c(p, resp), na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T))); abline(0,1)
sqrt((sum((resp-p)^2, na.rm=T))/length(p))

resp<-data$field_lai[data$train.1==0]
p <-predict(clf, data.matrix(data[data$train.1==0, feature.names]))
rmse<-sqrt((sum((resp-p)^2, na.rm=T))/length(p))
print(rmse)
plot(resp, p, ylim=c(0,max(c(p, resp),na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T)), pch=19,
     col=trans.color[c(8,1,2,4)][factor(data$site[data$train.1==0&data$site!="xtidb"])],
     cex.lab=1.5, cex.axis=1.2)
abline(0,1)
mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.1f",round(rmse,1))) ~ g ~ m^-2), cex=1.2)
#legend("bottomright",inset=c(0,0.1),legend=c("train", "test"), pch=c(19,1), bty="n", title=expression(underline(data)), title.adj = 0.3)
legend("bottomright", legend=unique(data$site[!is.na(data$site)&data$site!="xtidb"]), 
       title=expression(underline(site)),title.adj = 0.25,
       col=trans.color[c(8,1,2,4)][factor(unique(data$site[!is.na(data$site)&data$site!="xtidb"]))], pch=15, bty="n")

##aboveground biophysical results
#mae[n]
(sum((p-resp), na.rm=T))/length(resp[!is.na(resp)])
##RMSE
sqrt((sum((resp-p)^2, na.rm=T))/length(resp))
#nRMSE
(sqrt((sum((resp-p)^2, na.rm=T))/length(resp)))/(max(resp, na.rm=T)-min(resp, na.rm=T))
#COV RMSE (alternate nRMSE)
(sqrt((sum((resp-p)^2, na.rm=T))/length(resp)))/(mean(resp, na.rm=T))
#cor
cor(resp, p, use="pairwise.complete.obs")
#Ntest
nrow(data[!is.na(data$field_lai)&data$train.1==0,])
#Ntrain
nrow(data[!is.na(data$field_lai)&data$train.1==1,])



laimod1<-clf 
#clf<-laimod1

data<-data.frame(data)
data$predlai.l8<-(predict(laimod1, data.matrix(data[,feature.names])))
maxx<-max(c(data$field_lai, data$predlai.l8), na.rm=T)
minn<-min(c(data$field_lai, data$predlai.l8), na.rm=T)
sqrt((sum((data$field_lai-data$predlai.l8)^2, na.rm=T))/length(data$predlai.l8[!is.na(data$predlai.l8)]))

if(save_model==TRUE){
  ##save a copy of the model for future use
  saveRDS(laimod1, "output/xgb_lai_model.rds")
  xgb.save(laimod1,  "output/xgb_lai_model")
  laimod1<-xgb.load("output/xgb_lai_model")
  
  ##create a graph of measured vs predicted for training and testing data
  jpeg(file="results/xgb_l8_lai.jpg", res=200, units="in",
       height=5, width=5)
  par(mar=c(4,5,0.5,0.75), oma=c(0,0.2,0.5,0.5), mfrow=c(1,1))
  plot(data$field_lai, data$predlai.l8, 
       xlab="measured LAI", 
       ylab ="predicted LAI",
       ylim=c(0,maxx),xlim=c(0,maxx),  pch=c(1,19)[factor(data$train)], 
       cex=c(1,0.8)[factor(data$train)], col=trans.color[c(8,1,2,4)][factor(data$site)],
       cex.lab=1.5, cex.axis=1.2)
  abline(0,1)
  mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.1f",round(rmse,1)))), cex=1.2)
  legend("bottomright",inset=c(0,0.07),legend=c("train", "test"), pch=c(19,1), bty="n", title=expression(underline(data)), title.adj = 0.3)
  legend("bottomright", inset=c(0.2,0),legend=unique(data$site[!is.na(data$site)&data$site!="xtidb"]), 
         title=expression(underline(site)),title.adj = 0.25,
         col=trans.color[c(8,1,2,4)][factor(unique(data$site[!is.na(data$site)&data$site!="xtidb"]))], pch=15, bty="n")
  dev.off()
}


####CHL#####
tr<-data
tr<-dplyr::filter(tr, agm2.allom>10)

feature.names<-c("b1","b2", "b3","b4","b5","b6","b7", "gari", "ndmi", "pheno", "vari", "ndvi")
resp<-"chlorophyll"
tr<-data[data$train.1==1,c(resp, feature.names)]
tr<-tr[complete.cases(tr),]

traintask = TaskRegr$new(id = "id", backend = tr[,c(resp, feature.names)], target = resp)
traintask$select(feature.names)

##set up learner
learner = lrn("regr.xgboost", predict_type = "response",  booster="gbtree",  objective   = "reg:squarederror", nrounds=500)
#resamp=rsmp("holdout", ratio=0.85)
resamp=rsmp("cv", folds=inner)
measure = msr("regr.mse")
tuner = mlr3tuning::tnr("grid_search", resolution = 25)
terminator = trm("stagnation", iters=4, threshold=0.1)

instance = mlr3tuning::TuningInstanceSingleCrit$new(
  task = traintask,
  learner = learner,
  resampling = resamp,
  measure = measure,
  search_space = tune_ps,
  terminator = terminator
)
instance

##run the optimization
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
#instance$result_learner_param_vals

##store the best parameters back in the learner
learner$param_set$values = instance$result_learner_param_vals
learner$train(traintask)

##see the result, eg, the best configuration of parameters found during the iterations
param <- instance$result_x_domain

resp<-tr$chlorophyll
clf <- xgboost(data        = data.matrix(tr[,c(feature.names)]),
               label       = resp,
               booster="gbtree",
               nrounds     = 100,
               params = param,
               verbose=0)

mat <- xgb.importance (feature_names = feature.names,model = clf)
xgb.plot.importance (importance_matrix = mat[1:18]) 
p <-predict(clf, data.matrix(tr[, feature.names]))
plot(resp, p, ylim=c(0,max(c(p, resp), na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T))); abline(0,1)
sqrt((sum((resp-p)^2, na.rm=T))/length(p))

resp<-data$chlorophyll[data$train.1==0]
p <-predict(clf, data.matrix(data[data$train.1==0, feature.names]))
rmse<-sqrt((sum((resp-p)^2, na.rm=T))/length(p))
print(rmse)
plot(resp, p, ylim=c(0,max(c(p, resp),na.rm=T)),
     xlim=c(0,max(c(p, resp), na.rm=T)), pch=19,
     col=trans.color[c(8,1,2,4)][factor(data$site[data$train.1==0&data$site!="xtidb"])],
     cex.lab=1.5, cex.axis=1.2)
abline(0,1)
mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.2f",round(rmse,2))) ), cex=1.2)
#legend("bottomright",inset=c(0,0.1),legend=c("train", "test"), pch=c(19,1), bty="n", title=expression(underline(data)), title.adj = 0.3)
legend("bottomright", legend=unique(data$site[!is.na(data$site)&data$site!="xtidb"]), 
       title=expression(underline(site)),title.adj = 0.25,
       col=trans.color[c(8,1,2,4)][factor(unique(data$site[!is.na(data$site)&data$site!="xtidb"]))], pch=15, bty="n")

##aboveground biophysical results
#mae[n]
(sum((p-resp), na.rm=T))/length(resp[!is.na(resp)])
##RMSE
sqrt((sum((resp-p)^2, na.rm=T))/length(resp))
#nRMSE
(sqrt((sum((resp-p)^2, na.rm=T))/length(resp)))/(max(resp, na.rm=T)-min(resp, na.rm=T))
#COV RMSE (alternate nRMSE)
(sqrt((sum((resp-p)^2, na.rm=T))/length(resp)))/(mean(resp, na.rm=T))
#cor
cor(resp, p, use="pairwise.complete.obs")
#Ntest
nrow(data[!is.na(data$chlorophyll)&data$train.1==0,])
#Ntrain
nrow(data[!is.na(data$chlorophyll)&data$train.1==1,])

chlmod1<-clf 
#clf<-chlmod1

data<-data.frame(data)
data$predchl.l8<-(predict(chlmod1, data.matrix(data[,feature.names])))
maxx<-max(c(data$chlorophyll, data$predchl.l8), na.rm=T)
minn<-min(c(data$chlorophyll, data$predchl.l8), na.rm=T)
sqrt((sum((data$chlorophyll-data$predchl.l8)^2, na.rm=T))/length(data$predchl.l8[!is.na(data$predchl.l8)]))


if(save_model==TRUE){
  ##save a copy of the model for future use
  
  saveRDS(chlmod1, "output/xgb_chl_model.rds")
  xgb.save(chlmod1,  "output/xgb_chl_model")
  chlmod1<-xgb.load("output/xgb_chl_model")
  
  ##create a graph of measured vs predicted for training and testing data
  jpeg(file="results/xgb_l8_chl.jpg", res=200, units="in",
       height=5, width=5)
  par(mar=c(4,5,0.5,0.75), oma=c(0,0.2,0.5,0.5), mfrow=c(1,1))
  plot(data$chlorophyll, data$predchl.l8, 
       xlab =expression(paste("measured CHL mg ", g^-1, sep="")) , 
       ylab =expression(paste("predicted CHL mg ", g^-1, sep="")) , 
       ylim=c(0,maxx),xlim=c(0,maxx),  pch=c(1,19)[factor(data$train)], 
       cex=c(1,0.8)[factor(data$train)], col=trans.color[c(8,1,2,4)][factor(data$site)],
       cex.lab=1.5, cex.axis=1.2)
  abline(0,1)
  mtext(side=3, adj=0.05, line=-2, bquote(RMSE[test]: ~ .(sprintf("%.2f",round(rmse,2))) ~ mg ~ g^-1), cex=1.2)
  legend("bottomright",inset=c(0,0.07),legend=c("train", "test"), pch=c(19,1), bty="n", title=expression(underline(data)), title.adj = 0.3)
  legend("bottomright", inset=c(0.2,0),legend=unique(data$site[!is.na(data$site)&data$site!="xtidb"]), 
         title=expression(underline(site)),title.adj = 0.25,
         col=trans.color[c(8,1,2,4)][factor(unique(data$site[!is.na(data$site)&data$site!="xtidb"]))], pch=15, bty="n")
  dev.off()
}


##save a copy of the prediction results for later use
if(save_model==T){
  write_csv (data, path =  "output/xgb_predicted_biophysical_seagrant_landsat8.csv")
}

