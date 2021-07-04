##process newpixels to prepare for xgb bgbiomass model

library(R.utils); library(caret); library(tidyverse);library(RColorBrewer); library(mlr)
library(data.table);library(raster); library(xgboost); library(RcppRoll); library(lubridate)

##output and save new models and data sets?? Set to True/False below. 
##Note that True will overwrite the existing models and output datasets, so you'll want a backup if you don't intend this
save_model<-FALSE
save_model<-TRUE


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

##read in the greenup and soil temp data
newpixels.heat.l8<-read_csv("output/gce_pixels_landsat8_greenup_via_lst_and_air_lags.csv")
newpixels<-read_csv("output/xgb_predicted_biophysical_seagrant_landsat8_allpixels.csv")
newpixels$date<-as.Date(newpixels$date)
flood<-dplyr::select(newpixels, pix, flood_time)

newpixels.heat.l8$greenup<-as.Date(newpixels.heat.l8$greenup.soil)
newpixels.heat.l8$greenup<-as.Date(newpixels.heat.l8$greenup, origin="1970-01-01")

####
newpixels.daymet<-read_csv("output/daymet_monthly_climate_newpixels.csv")

###predict hybrid bg for newpixels
##smooth veg data series to estimate info on all valid landsat dates
dates<-seq.Date(as.Date("2013-04-15"), as.Date("2020-09-15"), by= "1 month")
newpixels<-newpixels[!is.na(newpixels$predn.l8),]

##predict once monthly estimate through approx, but much faster than for loop
##this allows dplyr to return a result that is greater > 1 per group
out<-newpixels %>% dplyr::group_by(pix, flood_time) %>% dplyr::filter(length(date)>25) %>%
  do(data.frame(date=dates,
                obs=length(.$date),
                year_obs=c(rep(length(.$date[.$year==2013]), length(dates[format(dates, "%Y")=="2013"])),
                           rep(length(.$date[.$year==2014]), length(dates[format(dates, "%Y")=="2014"])),
                           rep(length(.$date[.$year==2015]), length(dates[format(dates, "%Y")=="2015"])),
                           rep(length(.$date[.$year==2016]), length(dates[format(dates, "%Y")=="2016"])),
                           rep(length(.$date[.$year==2017]), length(dates[format(dates, "%Y")=="2017"])),
                           rep(length(.$date[.$year==2018]), length(dates[format(dates, "%Y")=="2018"])),
                           rep(length(.$date[.$year==2019]), length(dates[format(dates, "%Y")=="2019"])),
                           rep(length(.$date[.$year==2020]), length(dates[format(dates, "%Y")=="2020"]))),
                predn.l8= approx(as.numeric(.$date),  .$predn.l8, 
                                 xout=as.numeric(dates), method="linear")$y,
                predchl.l8= approx(as.numeric(.$date),  .$predchl.l8, 
                                   xout=as.numeric(dates), method="linear")$y,
                predlai.l8= approx(as.numeric(.$date),  .$predlai.l8, 
                                   xout=as.numeric(dates), method="linear")$y,
                predag.allom.l8= approx(as.numeric(.$date),  .$predag.allom.l8, 
                                        xout=as.numeric(dates), method="linear")$y,
                lst= approx(as.numeric(.$date),  .$lst, 
                            xout=as.numeric(dates), method="linear")$y,
                ndvi= approx(as.numeric(.$date),  .$ndvi, 
                            xout=as.numeric(dates), method="linear")$y
                
  )
  )
out$date<-as.Date(out$date, origin="1970-01-01")
out$year<-as.numeric(format(out$date, "%Y"))
out$mo<-as.numeric(format(out$date, "%m"))
out<-arrange(out, pix, year, mo)
newpixels<-out
newpixels.heat.l8<-ungroup(newpixels.heat.l8)
newpixels.heat.l8$year<-as.numeric(newpixels.heat.l8$year)
newpixels.heat.l8$key<-paste(newpixels.heat.l8$pix, newpixels.heat.l8$year)
newpixels.heat.l8<-dplyr::select(newpixels.heat.l8, -pix, -year)
newpixels$key<-paste(newpixels$pix, newpixels$year)
newpixels<-right_join(newpixels.heat.l8, newpixels, by="key")
newpixels<-dplyr::select(newpixels, -key)
newpixels<-arrange(newpixels, pix, year, mo)
newpixels$lpredn<-log(newpixels$predn.l8)
newpixels$lpredag<-log(newpixels$predag.allom.l8)
newpixels$lpredlai<-log(newpixels$predlai.l8)
newpixels$greendoy<-as.numeric(format(newpixels$greenup, "%j"))
newpixels$doy<-as.numeric(format(newpixels$date, "%j"))
newpixels$growingday<-newpixels$doy-newpixels$greendoy
newpixels$growingday<-ifelse(newpixels$doy>305|newpixels$doy<newpixels$greendoy, 0,newpixels$growingday)

newpixels<-newpixels %>% dplyr::group_by( pix) %>%
  dplyr::mutate(photo=predn.l8/100*predag.allom.l8,
                #photo=predchl.l8*predlai.l8,
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
newpixels$greenupmo<-as.numeric(format(round_date(newpixels$greenup, unit="month"), "%m"))
yrs<-unique(newpixels$year); pixels<-unique(newpixels$pix)
newpixels$startag<-NA
for (j in seq_along(yrs)){
  for (i in seq_along(pixels)){
    startag<-newpixels$predag.allom.l8[newpixels$year==yrs[j]&newpixels$pix==pixels[i]&newpixels$mo==newpixels$greenupmo] 
    newpixels$startag[newpixels$year==yrs[j]&newpixels$pix==pixels[i]]<-
      rep(startag, length(newpixels$predag.allom.l8[newpixels$year==yrs[j]&newpixels$pix==pixels[i]]))
    
    startn<-newpixels$predn.l8[newpixels$year==yrs[j]&newpixels$pix==pixels[i]&newpixels$mo==newpixels$greenupmo] 
    newpixels$startn[newpixels$year==yrs[j]&newpixels$pix==pixels[i]]<-
      rep(startn, length(newpixels$predn.l8[newpixels$year==yrs[j]&newpixels$pix==pixels[i]]))
    
    startlai<-newpixels$predlai.l8[newpixels$year==yrs[j]&newpixels$pix==pixels[i]&newpixels$mo==newpixels$greenupmo] 
    newpixels$startlai[newpixels$year==yrs[j]&newpixels$pix==pixels[i]]<-
      rep(startlai, length(newpixels$predlai.l8[newpixels$year==yrs[j]&newpixels$pix==pixels[i]]))
    
    startchl<-newpixels$predchl.l8[newpixels$year==yrs[j]&newpixels$pix==pixels[i]&newpixels$mo==newpixels$greenupmo] 
    newpixels$startchl[newpixels$year==yrs[j]&newpixels$pix==pixels[i]]<-
      rep(startchl, length(newpixels$predchl.l8[newpixels$year==yrs[j]&newpixels$pix==pixels[i]]))
    
    startphoto<-newpixels$photo[newpixels$year==yrs[j]&newpixels$pix==pixels[i]&newpixels$mo==newpixels$greenupmo] 
    newpixels$startphoto[newpixels$year==yrs[j]&newpixels$pix==pixels[i]]<-
      rep(startphoto, length(newpixels$photo[newpixels$year==yrs[j]&newpixels$pix==pixels[i]]))
    ##correct the early months to reflect the previous growing season
    if(yrs[j]>2013){
      newpixels$startag[newpixels$year==yrs[j]&newpixels$pix==pixels[i]&newpixels$doy<newpixels$greendoy]<-newpixels$predag.allom.l8[newpixels$year==(yrs[j]-1)&newpixels$pix==pixels[i]&newpixels$mo==newpixels$greenupmo] 
      newpixels$startn[newpixels$year==yrs[j]&newpixels$pix==pixels[i]&newpixels$doy<newpixels$greendoy]<-newpixels$predn.l8[newpixels$year==(yrs[j]-1)&newpixels$pix==pixels[i]&newpixels$mo==newpixels$greenupmo] 
      newpixels$startlai[newpixels$year==yrs[j]&newpixels$pix==pixels[i]&newpixels$doy<newpixels$greendoy]<-newpixels$predlai.l8[newpixels$year==(yrs[j]-1)&newpixels$pix==pixels[i]&newpixels$mo==newpixels$greenupmo] 
      newpixels$startchl[newpixels$year==yrs[j]&newpixels$pix==pixels[i]&newpixels$doy<newpixels$greendoy]<-newpixels$predchl.l8[newpixels$year==(yrs[j]-1)&newpixels$pix==pixels[i]&newpixels$mo==newpixels$greenupmo] 
      newpixels$startphoto[newpixels$year==yrs[j]&newpixels$pix==pixels[i]&newpixels$doy<newpixels$greendoy]<-newpixels$photo[newpixels$year==(yrs[j]-1)&newpixels$pix==pixels[i]&newpixels$mo==newpixels$greenupmo] 
    }
  }
}

newpixels$diff_ag_growing<- newpixels$predag.allom.l8-newpixels$startag
newpixels$diff_N_growing<- newpixels$predn.l8-newpixels$startn
newpixels$diff_chl_growing<- newpixels$predchl.l8-newpixels$startchl
newpixels$diff_lai_growing<- newpixels$predlai.l8-newpixels$startlai
newpixels$diff_photo_growing<- newpixels$photo-newpixels$startphoto
##create growing season sums: the sum of the proceeding season assigned to the next winter,
##later we fill in the non-winter months with the rolling mean of 4 months
tmp<- newpixels %>% group_by(pix, year) %>% filter(mo==10)%>%
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
newpixels$key<-paste(newpixels$mo, newpixels$year, newpixels$pix)
newpixels<-data.table(newpixels)
setkey(newpixels, key)
newpixels<-merge(newpixels, tmp, by="key", all.x=T)

newpixels<-data.frame(newpixels)
remove(tmp); remove(tmp2);remove(tmp3)

##fill in the growing season months with rolling sums, winter months are the sum of the previous season
newpixels$growag<-ifelse(is.na(newpixels$growag), newpixels$roll_ag_4, newpixels$growag)
newpixels$growN<-ifelse(is.na(newpixels$growN), newpixels$roll_N_4, newpixels$growN)
newpixels$growlai<-ifelse(is.na(newpixels$growlai), newpixels$roll_lai_4, newpixels$growlai)
newpixels$growchl<-ifelse(is.na(newpixels$growchl), newpixels$roll_chl_4, newpixels$growchl)
newpixels$growphoto<-ifelse(is.na(newpixels$growphoto), newpixels$roll_photo_4, newpixels$growphoto)

##extract dem elevation from dem ascii layer by utm coordinates
#library(raster)
x<-str_split(newpixels$pix, pattern=" ", simplify = T)
newpixels$utm_east<-as.numeric(x[,1])
newpixels$utm_north<-as.numeric(x[,2])

plot.locations<-data.frame(newpixels[,c("utm_east", "utm_north")])
dem<-raster("/home/jessica/UGA/Projects/seagrant/data/ModDEM_ascii/mod_dem.txt")

plot.dem<-raster::extract(dem,plot.locations)
newpixels$elevation<-plot.dem
newpixels$mo<-as.numeric(format(newpixels$date, "%m"))
newpixels<-newpixels[!is.na(newpixels$date),]

##merge with tide data
newpixels$key<-paste(newpixels$year, newpixels$mo)
newpixels<-merge(newpixels, tide, by="key", all.x=T, all.y=T)
newpixels$local_hiwater<-(newpixels$hiwater-newpixels$elevation)/(newpixels$hiwater-newpixels$water)*newpixels$flood_time#
newpixels$local_lowater<-(newpixels$elevation-newpixels$lowater)/(newpixels$water-newpixels$lowater)*(1-newpixels$flood_time)
newpixels<-arrange(newpixels, pix, date)
newpixels<- newpixels %>% group_by(pix) %>% mutate(
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
newpixels<-ungroup(newpixels)
newpixels$loc_elevation<-(newpixels$elevation-newpixels$water)/(newpixels$hiwater-newpixels$water)

##merge daymet data
newpixels$key<-paste(newpixels$date, newpixels$pix)
newpixels<-merge(newpixels, newpixels.daymet, by="key", all.x=T)

if(save_model==T){
  write_csv(newpixels, "output/xgb_processed_seagrant_landsat8_allpixels.csv")
}