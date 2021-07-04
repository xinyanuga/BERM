###Build xgboost belowground biomass model from predicted aboveground biophysical data at ground truth locations
##this takes the data created from xgboost_landsat_biophysical.r

library(R.utils); library(caret); library(tidyverse);library(RColorBrewer);
library(data.table);library(raster); library(xgboost); library(RcppRoll); library(lubridate)

setwd("/home/jessica/Documents/UGA/Projects/BERM")

##set up some colors for plots
color<-brewer.pal(8, "Dark2")
trans.color<-addalpha(color,0.9)


##read in estimate of percent of flooded obs per pixel, as estimated by the Landsat 8 super flood model (RF)
flood<-read_csv("output/flooded_percent_from_landsat_seagrant_pixels.csv")

##read in the estimated green up day of year, lst and soil temp data
heat.l8<-read_csv("output/seagrant_pixels_landsat8_heat_processed.csv")

####read in gridded climate data from Daymet, aquired via GEE and processed to monthly estimates by the process_daymet.r script
daymet<-read_csv("output/daymet_monthly_climate_seagrant_plots.csv")

##read in the predictions for the aboveground biophysical models
data<-read_csv ("output/xgb_predicted_biophysical_seagrant_landsat8.csv")
data$date<-as.Date(data$date)
data<-data[data$site!="xtidb",]
heat.l8$greenup<-as.Date(heat.l8$greenup.soil)
heat.l8$greenup<-as.Date(heat.l8$greenup, origin="1970-01-01")


####Ground-truth data for field observations####
##Read in location data for vegetation ground truth plots
plots<-read.csv("data/biomass_plots.csv", header=T)
##not all dates have all pixels, so we need to make the averages here in order to avoid duplicate elevations later
pixels<- plots %>% group_by(pix) %>% summarise(elevation=mean(elevation), 
                                               utm_east=mean(utm_east), 
                                               utm_north=mean(utm_north), 
                                               lat=mean(lat), long=mean(long))
plots<-dplyr::select(plots, plot, pix, elevation)


##Creating middle of month approximations for all data, so that we have an even data set with one monthly observation
#i<-1; j<-1; k<-1
var<-c("percentN", "bgm2.areascaled","agm2.core.areascaled", "bgm2.core.stemscaled",
       "agm2.core.stemscaled", "agm2.allom", "bgm2.allometric.rootshoot",
       "bgm2.allometric.rootgreenshoot", "agm2.allom.new", "bgm2.allo.new.rootshoot",
       "bgm2.allo.new.rootgreenshoot", "chlorophyll", "field_lai", 
       "predag.allom.l8", "predn.l8", "predchl.l8", "predlai.l8", "lst", "ndvi")
loc<-unique(data$pix)
dates<-seq.Date(as.Date("2013-04-15"), as.Date("2020-10-01"), by= "1 month")
out<-data.frame(matrix(ncol=length(var)+2, nrow=5000))
names(out)<-c("date", "pix", var)
loop<-1
data<-data.frame(data)
for (j in seq_along(loc)){
  for (i in seq_along(var)){
    x<-data[data$pix==loc[j],]
    x<-x[,c("date", var[i])]
    x<-x[complete.cases(x),]
    if(nrow(x)>1)   {
      smoothed <- approx(as.numeric(x$date),x[,2], xout=as.numeric(dates),
                         method="linear"    )
      num<-length(smoothed$x)
      out[loop:(loop+num-1),"pix"]<-rep(loc[j], num)
      out[loop:(loop+num-1),"date"]<-as.Date(smoothed$x, origin = "1970-01-01")
      out[loop:(loop+num-1), var[i]]<-(smoothed$y)
    }
  }
  loop<-loop+num    
}

out<-out[!is.na(out$date),]
out$date<-as.Date(out$date, origin="1970-01-01")
out$year<-as.numeric(format(out$date, "%Y"))
out$mo<-as.numeric(format(out$date, "%m"))
out$site<-substr(out$pix, 1,2)
out$site[!(out$site %in% c("fb", "sk", "ug"))]<- "fa"
out$site<-ifelse(out$site=="fb", "fluxb", out$site)
out$site<-ifelse(out$site=="fa", "fluxa", out$site)
out$site<-ifelse(out$site=="ug", "ugami", out$site)
out$site<-ifelse(out$site=="sk", "skida", out$site)
out<-arrange(out, pix, year, mo)


##combine with greenup info
heat.l8<-ungroup(heat.l8)
heat.l8$key<-paste(heat.l8$pix, heat.l8$year)
heat.l8<-dplyr::select(heat.l8, -pix, -year)
out$key<-paste(out$pix, out$year)
out<-full_join(heat.l8, out, by="key", all.y=T)
out<-dplyr::select(out, -key)

##grab elevation data to add to the new dataset
elevation<-unique(dplyr::select(data, pix, elevation))

##now we don't need the old data set, so we can write over it
data<-left_join(out, elevation)
data<-arrange(data, pix, year, mo)
data$lpredn<-log(data$predn.l8)
data$lpredag<-log(data$predag.allom.l8)
data$lpredlai<-log(data$predlai.l8)
data$greendoy<-as.numeric(format(data$greenup, "%j"))
data$doy<-as.numeric(format(data$date, "%j"))
##create growing day variable= days since green up during the growing season, 0 during winter
data$growingday<-data$doy-data$greendoy
data$growingday<-ifelse(data$doy>305|data$doy<data$greendoy, 0,data$growingday)

##prepare potential predictor variables for the Extreme Gradient Boosting model
##we want combinations of antecedent conditions represented in the prediction set
data<-data %>% dplyr::group_by( pix) %>%
  dplyr::mutate(photo=predn.l8/100*predag.allom.l8,
    lagg_photo_1=lag(photo),
    diff_photo_1=c(NA, diff(photo, lag=1)),
    diff_photo_2=c(rep(NA, 2), diff(photo, lag=2)),
    diff_photo_3=c(rep(NA, 3),diff(photo, lag=3)),
    diff_photo_4=c(rep(NA, 4), diff(photo, lag=4)),
    diff_photo_5=c(rep(NA, 5),diff(photo, lag=5)),
    diff_photo_6=c(rep(NA, 6), diff(photo, lag=6)),
    roll_photo_2=RcppRoll::roll_meanr(photo, n=2, fill=NA, na.rm=T),
    roll_photo_3=RcppRoll::roll_meanr(photo, n=3, fill=NA, na.rm=T),
    roll_photo_4=RcppRoll::roll_meanr(photo, n=4, fill=NA, na.rm=T),
    roll_photo_5=RcppRoll::roll_meanr(photo, n=5, fill=NA, na.rm=T),
    roll_photo_6=RcppRoll::roll_meanr(photo, n=6, fill=NA, na.rm=T),
    lagg_lai_1=lag(predlai.l8),
    diff_lai_1=c(NA, diff(predlai.l8, lag=1)),
    diff_lai_2=c(rep(NA, 2), diff(predlai.l8, lag=2)),
    diff_lai_3=c(rep(NA, 3),diff(predlai.l8, lag=3)),
    diff_lai_4=c(rep(NA, 4), diff(predlai.l8, lag=4)),
    diff_lai_5=c(rep(NA, 5),diff(predlai.l8, lag=5)),
    diff_lai_6=c(rep(NA, 6), diff(predlai.l8, lag=6)),
    roll_lai_2=RcppRoll::roll_meanr(predlai.l8, n=2, fill=NA, na.rm=T),
    roll_lai_3=RcppRoll::roll_meanr(predlai.l8, n=3, fill=NA, na.rm=T),
    roll_lai_4=RcppRoll::roll_meanr(predlai.l8, n=4, fill=NA, na.rm=T),
    roll_lai_5=RcppRoll::roll_meanr(predlai.l8, n=5, fill=NA, na.rm=T),
    roll_lai_6=RcppRoll::roll_meanr(predlai.l8, n=6, fill=NA, na.rm=T),
    lagg_ag_1=lag(predag.allom.l8),
    diff_ag_1=c(NA, diff(predag.allom.l8, lag=1)),
    diff_ag_2=c(rep(NA, 2), diff(predag.allom.l8, lag=2)),
    diff_ag_3=c(rep(NA, 3),diff(predag.allom.l8, lag=3)),
    diff_ag_4=c(rep(NA, 4), diff(predag.allom.l8, lag=4)),
    diff_ag_5=c(rep(NA, 5),diff(predag.allom.l8, lag=5)),
    diff_ag_6=c(rep(NA, 6), diff(predag.allom.l8, lag=6)),
    roll_ag_2=RcppRoll::roll_meanr(predag.allom.l8, n=2, fill=NA, na.rm=T),
    roll_ag_3=RcppRoll::roll_meanr(predag.allom.l8, n=3, fill=NA, na.rm=T),
    roll_ag_4=RcppRoll::roll_meanr(predag.allom.l8, n=4, fill=NA, na.rm=T),
    roll_ag_5=RcppRoll::roll_meanr(predag.allom.l8, n=5, fill=NA, na.rm=T),
    roll_ag_6=RcppRoll::roll_meanr(predag.allom.l8, n=6, fill=NA, na.rm=T),
    lagg_lst_1=lag(lst),
    diff_lst_1=c(NA, diff(lst, lag=1)),
    diff_lst_2=c(rep(NA, 2), diff(lst, lag=2)),
    diff_lst_3=c(rep(NA, 3),diff(lst, lag=3)),
    diff_lst_4=c(rep(NA, 4), diff(lst, lag=4)),
    diff_lst_5=c(rep(NA, 5),diff(lst, lag=5)),
    diff_lst_6=c(rep(NA, 6), diff(lst, lag=6)),
    roll_lst_2=RcppRoll::roll_meanr(lst, n=2, fill=NA, na.rm=T),
    roll_lst_3=RcppRoll::roll_meanr(lst, n=3, fill=NA, na.rm=T),
    roll_lst_4=RcppRoll::roll_meanr(lst, n=4, fill=NA, na.rm=T),
    roll_lst_5=RcppRoll::roll_meanr(lst, n=5, fill=NA, na.rm=T),
    roll_lst_6=RcppRoll::roll_meanr(lst, n=6, fill=NA, na.rm=T),
    lagg_N_1=lag(predn.l8),
    diff_N_1=c(NA, diff(predn.l8, lag=1)),
    diff_N_2=c(rep(NA, 2), diff(predn.l8, lag=2)),
    diff_N_3=c(rep(NA, 3),diff(predn.l8, lag=3)),
    diff_N_4=c(rep(NA, 4), diff(predn.l8, lag=4)),
    diff_N_5=c(rep(NA, 5),diff(predn.l8, lag=5)),
    diff_N_6=c(rep(NA, 6), diff(predn.l8, lag=6)),
    roll_N_2=RcppRoll::roll_meanr(predn.l8, n=2, fill=NA, na.rm=T),
    roll_N_3=RcppRoll::roll_meanr(predn.l8, n=3, fill=NA, na.rm=T),
    roll_N_4=RcppRoll::roll_meanr(predn.l8, n=4, fill=NA, na.rm=T),
    roll_N_5=RcppRoll::roll_meanr(predn.l8, n=5, fill=NA, na.rm=T),
    roll_N_6=RcppRoll::roll_meanr(predn.l8, n=6, fill=NA, na.rm=T),
    lagg_chl_1=lag(predchl.l8),
    diff_chl_1=c(NA, diff(predchl.l8, lag=1)),
    diff_chl_2=c(rep(NA, 2), diff(predchl.l8, lag=2)),
    diff_chl_3=c(rep(NA, 3),diff(predchl.l8, lag=3)),
    diff_chl_4=c(rep(NA, 4), diff(predchl.l8, lag=4)),
    diff_chl_5=c(rep(NA, 5),diff(predchl.l8, lag=5)),
    diff_chl_6=c(rep(NA, 6), diff(predchl.l8, lag=6)),
    roll_chl_2=RcppRoll::roll_meanr(predchl.l8, n=2, fill=NA, na.rm=T),
    roll_chl_3=RcppRoll::roll_meanr(predchl.l8, n=3, fill=NA, na.rm=T),
    roll_chl_4=RcppRoll::roll_meanr(predchl.l8, n=4, fill=NA, na.rm=T),
    roll_chl_5=RcppRoll::roll_meanr(predchl.l8, n=5, fill=NA, na.rm=T),
    roll_chl_6=RcppRoll::roll_meanr(predchl.l8, n=6, fill=NA, na.rm=T), 
    lagg_ndvi_1=lag(ndvi),
    diff_ndvi_1=c(NA, diff(ndvi, lag=1)),
    diff_ndvi_2=c(rep(NA, 2), diff(ndvi, lag=2)),
    diff_ndvi_3=c(rep(NA, 3),diff(ndvi, lag=3)),
    diff_ndvi_4=c(rep(NA, 4), diff(ndvi, lag=4)),
    diff_ndvi_5=c(rep(NA, 5),diff(ndvi, lag=5)),
    diff_ndvi_6=c(rep(NA, 6), diff(ndvi, lag=6)),
    roll_ndvi_2=RcppRoll::roll_meanr(ndvi, n=2, fill=NA, na.rm=T),
    roll_ndvi_3=RcppRoll::roll_meanr(ndvi, n=3, fill=NA, na.rm=T),
    roll_ndvi_4=RcppRoll::roll_meanr(ndvi, n=4, fill=NA, na.rm=T),
    roll_ndvi_5=RcppRoll::roll_meanr(ndvi, n=5, fill=NA, na.rm=T),
    roll_ndvi_6=RcppRoll::roll_meanr(ndvi, n=6, fill=NA, na.rm=T)
  )

##calculate differences from time x since greenup, 
##where start[var] is the data value of var during the month of greenup for that pixel/year
data$greenupmo<-as.numeric(format(round_date(data$greenup, unit="month"), "%m"))
yrs<-unique(data$year[data$year>2013]); pixels<-unique(data$pix)
data$startag<-NA
#j<-1; i<-1
for (j in seq_along(yrs)){
  for (i in seq_along(pixels)){
    startag<-data$predag.allom.l8[data$year==yrs[j]&data$pix==pixels[i]&data$mo==data$greenupmo] 
    data$startag[data$year==yrs[j]&data$pix==pixels[i]]<-
      rep(startag, length(data$predag.allom.l8[data$year==yrs[j]&data$pix==pixels[i]]))
    
    startn<-data$predn.l8[data$year==yrs[j]&data$pix==pixels[i]&data$mo==data$greenupmo] 
    data$startn[data$year==yrs[j]&data$pix==pixels[i]]<-
      rep(startn, length(data$predn.l8[data$year==yrs[j]&data$pix==pixels[i]]))
    
    startlai<-data$predlai.l8[data$year==yrs[j]&data$pix==pixels[i]&data$mo==data$greenupmo] 
    data$startlai[data$year==yrs[j]&data$pix==pixels[i]]<-
      rep(startlai, length(data$predlai.l8[data$year==yrs[j]&data$pix==pixels[i]]))
    
    startchl<-data$predchl.l8[data$year==yrs[j]&data$pix==pixels[i]&data$mo==data$greenupmo] 
    data$startchl[data$year==yrs[j]&data$pix==pixels[i]]<-
      rep(startchl, length(data$predchl.l8[data$year==yrs[j]&data$pix==pixels[i]]))
    
    startphoto<-data$photo[data$year==yrs[j]&data$pix==pixels[i]&data$mo==data$greenupmo] 
    data$startphoto[data$year==yrs[j]&data$pix==pixels[i]]<-
      rep(startphoto, length(data$photo[data$year==yrs[j]&data$pix==pixels[i]]))
    
    ##correct the early months before green up to reflect values from the previous growing season
        if(yrs[j]>2013){
          data$startag[data$year==yrs[j]&data$pix==pixels[i]&data$doy<data$greendoy]<-data$predag.allom.l8[data$year==(yrs[j]-1)&data$pix==pixels[i]&data$mo==data$greenupmo] 
          data$startn[data$year==yrs[j]&data$pix==pixels[i]&data$doy<data$greendoy]<-data$predn.l8[data$year==(yrs[j]-1)&data$pix==pixels[i]&data$mo==data$greenupmo] 
          data$startlai[data$year==yrs[j]&data$pix==pixels[i]&data$doy<data$greendoy]<-data$predlai.l8[data$year==(yrs[j]-1)&data$pix==pixels[i]&data$mo==data$greenupmo] 
          data$startchl[data$year==yrs[j]&data$pix==pixels[i]&data$doy<data$greendoy]<-data$predchl.l8[data$year==(yrs[j]-1)&data$pix==pixels[i]&data$mo==data$greenupmo] 
          data$startphoto[data$year==yrs[j]&data$pix==pixels[i]&data$doy<data$greendoy]<-data$photo[data$year==(yrs[j]-1)&data$pix==pixels[i]&data$mo==data$greenupmo] 
    }
  }
}

data$diff_ag_growing<- data$predag.allom.l8-data$startag
data$diff_N_growing<- data$predn.l8-data$startn
data$diff_chl_growing<- data$predchl.l8-data$startchl
data$diff_lai_growing<- data$predlai.l8-data$startlai
data$diff_photo_growing<- data$photo-data$startphoto

##create growing season sums: the sum of the proceeding season assigned to the next winter,
##later we fill in the non-winter months with the rolling mean of 4 months
tmp<- data %>% group_by(pix, year) %>% filter(mo==10)%>%
  summarise(growag= roll_ag_2,
            growN=roll_N_2,
            growchl= roll_chl_2,
            growlai= roll_lai_2,
            growphoto=roll_photo_2
  )

tmp$key<-paste("10", tmp$year, tmp$pix)
tmp2<-tmp
tmp2$key<-paste("11", tmp$year, tmp$pix)
tmp3<-tmp
tmp3$key<-paste("12", tmp$year, tmp$pix)
tmp2<-rbind(tmp2,tmp3)
tmp3<-tmp
tmp3$key<-paste("1", tmp$year-1, tmp$pix)
tmp2<-rbind(tmp2,tmp3)
tmp3<-tmp
tmp3$key<-paste("2", tmp$year-1, tmp$pix)
tmp2<-rbind(tmp2,tmp3)
tmp3<-tmp
tmp3$key<-paste("3", tmp$year-1, tmp$pix)
tmp2<-rbind(tmp2,tmp3)
tmp3<-tmp
tmp3$key<-paste("4", tmp$year-1, tmp$pix)
tmp2<-rbind(tmp2,tmp3)
tmp<-rbind(tmp,tmp2)
tmp<-ungroup(tmp)
tmp<-dplyr::select(tmp, -pix, -year)
tmp<-data.table(tmp)
setkey(tmp, key)
data$key<-paste(data$mo, data$year, data$pix)
data<-data.table(data)
setkey(data, key)
data<-merge(data, tmp, by="key", all.x=T)

data<-data.frame(data)
remove(tmp); remove(tmp2);remove(tmp3)

###correct previous growing season var (a variable important in winter) to be the 
##rolling mean of the last 4 months if it isn't winter, thus it always reflects a long-term trend
data$growag<-ifelse(is.na(data$growag), data$roll_ag_4, data$growag)
data$growN<-ifelse(is.na(data$growN), data$roll_N_4, data$growN)
data$growlai<-ifelse(is.na(data$growlai), data$roll_lai_4, data$growlai)
data$growchl<-ifelse(is.na(data$growchl), data$roll_chl_4, data$growchl)
data$growphoto<-ifelse(is.na(data$growphoto), data$roll_photo_4, data$growphoto)


###Merge in Daymet monthly gridded climate data
data$key<-paste(data$date, data$pix)
data<-merge(data, daymet, by="key", all.x=T)

data<-merge(data, flood, by="pix", all.x=T)

###Monthly benchmarked water level means from fort pulaski, not detrended data
tide<-read.csv("data/fort_pulaski_monthly_mean_sea_level_NAVD_CO-OPS_8670870_wl_2013_2020.csv", stringsAsFactors = F)
##name the vars something nicer
tide$date<-as.Date(tide$Date, format="%Y/%m/%d")
tide$maxwater<-tide$Highest
tide$minwater<-tide$'Lowest..m.'
tide$hiwater<-tide$MHHW..m.
tide$lowater<-tide$MLLW..m.
tide$water<-tide$'MSL..m.'
tide$year<-as.numeric(format(tide$date, "%Y"))
tide$mo<-as.numeric(format(tide$date, "%m"))
tide<-dplyr::select(tide, date, year, mo, water, maxwater, minwater, hiwater, lowater)
tide<-tide[tide$year>2012,]
tide$key<-paste(tide$year, tide$mo)
tide<-dplyr::select(tide, key, water, hiwater, maxwater,lowater)

data$key<-paste(data$year, data$mo)
data<-merge(data, tide, by="key", all.x=T, all.y=T)

##create approx pixel flood depth intensity by subtracting pixel elevation from Mean High High Water 
data$local_hiwater<-(data$hiwater-data$elevation)/(data$hiwater-data$water)*data$flood_time
##create approx pixel lack of tidal flushing intensity by subtracting Mean Low Low  Water from pixel elevation 
data$local_lowater<-(data$elevation-data$lowater)/(data$water-data$lowater)*(1-data$flood_time)

data<-unique(data)
data<-ungroup(data)

##process antecedant conditions for tide/water level data too
data<-arrange(data, pix, date)
data<- data %>% group_by(pix) %>% mutate(
  lagg_lochi_1=lag(local_hiwater),
  roll_lochi_0=(local_hiwater),
  roll_lochi_1=RcppRoll::roll_sumr(lagg_lochi_1, n=2, fill=NA),
  roll_lochi_2=RcppRoll::roll_sumr(local_hiwater, n=2, fill=NA),
  roll_lochi_3=RcppRoll::roll_sumr(local_hiwater, n=3, fill=NA),
  roll_lochi_4=RcppRoll::roll_sumr(local_hiwater, n=4, fill=NA),
  roll_lochi_5=RcppRoll::roll_sumr(local_hiwater, n=5, fill=NA),
  lagg_loclo_1=lag(local_lowater),
  roll_loclo_0=(local_lowater),
  roll_loclo_1=RcppRoll::roll_sumr(lagg_loclo_1, n=2, fill=NA),
  roll_loclo_2=RcppRoll::roll_sumr(local_lowater, n=2, fill=NA),
  roll_loclo_3=RcppRoll::roll_sumr(local_lowater, n=3, fill=NA),
  roll_loclo_4=RcppRoll::roll_sumr(local_lowater, n=4, fill=NA),
  roll_loclo_5=RcppRoll::roll_sumr(local_lowater, n=5, fill=NA),
)
data<-ungroup(data)
data$id<-row.names(data)
data$year<-as.numeric(format(data$date, "%Y"))
data$mo<-(format(data$date, "%m"))
data$temp_diff<-(data$tmax_mean+data$tmin_mean)/2-data$lst
data$loc_elevation<-(data$elevation-data$water)/(data$hiwater-data$water)

##now we're read to actually build the model
##proceed to next script, model_xgb_belowground_biomass.r
