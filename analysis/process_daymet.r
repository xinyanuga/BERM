###xgboost belowground biomass
library(tidyverse);

setwd("/home/jessica/Documents/UGA/Projects/BERM")

##process daymet data
##make date equivalent to estimated veg data dates


##Read in the daymet data acquired from google earth engine for the site level landsat pixel locations
daymet<-read_csv("data/daymet_daily_climate_seagrantpts_2013_2020.csv")

##process date and location variables
daymet$date<-as.Date(daymet$date)
daymet$mo<-format(daymet$date, "%m")
daymet$yr<-format(daymet$date, "%Y")
daymet$site<-substr(daymet$plot, 1,2)
daymet$site[!(daymet$site %in% c("fb", "sk", "ug"))]<- "fa"
##remove variables we don't need and create estimates summarized by month and plot
daymet<-dplyr::select(daymet, -'system:index', -lat, -long, -utm_east, -utm_north, -rtk_elev, -dem_rtk_el, -zone, -.geo, -date, -swe)
daymet<-daymet %>% group_by(yr, mo, plot, site) %>% summarise_all(funs(mean, sum))
##assign the monthly estimate to the middle of the month
daymet$date<-as.Date(paste(daymet$yr, daymet$mo, "15", sep="-"))

##name the pixels so that you can group by them
daymet$pix<-"xx"
daymet$pix[(daymet$plot %in% c("sk1", "sk2", "sk3"))]<-"ska"
daymet$pix<-ifelse(daymet$plot %in% c("sk4", "sk5", "sk6"),"skb",daymet$pix)
daymet$pix<-ifelse(daymet$plot %in% c("sk7", "sk8", "sk9"),"skc",daymet$pix)
daymet$pix<-ifelse(daymet$plot %in% c("fb1", "fb2", "fb3"),"fba",daymet$pix)
daymet$pix<-ifelse(daymet$plot %in% c("fb4", "fb5", "fb6"),"fbb",daymet$pix)
daymet$pix<-ifelse(daymet$plot %in% c("fb7", "fb8", "fb9"),"fbc",daymet$pix)
daymet$pix<-ifelse(daymet$plot %in% c("ug1", "ug2", "ug3"),"uga",daymet$pix)
daymet$pix<-ifelse(daymet$plot %in% c("ug4", "ug5", "ug6"),"ugb",daymet$pix)
daymet$pix<-ifelse(daymet$plot %in% c("ug7", "ug8", "ug9"),"ugc",daymet$pix)
daymet$pix<-ifelse(daymet$plot %in% c("t7","t8"),"faf",daymet$pix)
daymet$pix<-ifelse(daymet$plot %in% c("t1", "t2", "t3", "t4"),"faa",daymet$pix)
daymet$pix<-ifelse(daymet$plot %in% c("t5", "t6", "m1", "m2","m3"),"fab",daymet$pix)
daymet$pix<-ifelse(daymet$plot %in% c("m4", "m5", "m6"),"fac",daymet$pix)
daymet$pix<-ifelse(daymet$plot %in% c("s1", "s2", "s6"),"fad",daymet$pix)
daymet$pix<-ifelse(daymet$plot %in% c("s3", "s4", "s5"),"fae",daymet$pix)
daymet<-ungroup(daymet)
daymet<-dplyr::select(daymet, -yr, -mo, -plot,-site)
daymet<-daymet %>% group_by(date, pix) %>% summarise_all(funs(mean))

##calculate PAR from solar radiation via standard formula
daymet$par_sum<-2.114*daymet$srad_sum
daymet$par_tot<-daymet$dayl_sum*daymet$par_sum/1000000

##create derived varialbes that we might use in the Machine Learning Model
daymet<-arrange(daymet, pix, date)
daymet<-daymet %>% dplyr::group_by( pix) %>%
  dplyr::mutate(lagg_tmax_1=lag(tmax_mean),
                diff_tmax_1=c(NA, diff(tmax_mean, lag=1)),
                diff_tmax_2=c(rep(NA, 2), diff(tmax_mean, lag=2)),
                diff_tmax_3=c(rep(NA, 3),diff(tmax_mean, lag=3)),
                diff_tmax_4=c(rep(NA, 4), diff(tmax_mean, lag=4)),
                diff_tmax_5=c(rep(NA, 5),diff(tmax_mean, lag=5)),
                diff_tmax_6=c(rep(NA, 6), diff(tmax_mean, lag=6)),
                roll_tmax_2=RcppRoll::roll_meanr(tmax_mean, n=2, fill=NA, na.rm=T),
                roll_tmax_3=RcppRoll::roll_meanr(tmax_mean, n=3, fill=NA, na.rm=T),
                roll_tmax_4=RcppRoll::roll_meanr(tmax_mean, n=4, fill=NA, na.rm=T),
                roll_tmax_5=RcppRoll::roll_meanr(tmax_mean, n=5, fill=NA, na.rm=T),
                roll_tmax_6=RcppRoll::roll_meanr(tmax_mean, n=6, fill=NA, na.rm=T),
                lagg_tmin_1=lag(tmin_mean),
                diff_tmin_1=c(NA, diff(tmin_mean, lag=1)),
                diff_tmin_2=c(rep(NA, 2), diff(tmin_mean, lag=2)),
                diff_tmin_3=c(rep(NA, 3),diff(tmin_mean, lag=3)),
                diff_tmin_4=c(rep(NA, 4), diff(tmin_mean, lag=4)),
                diff_tmin_5=c(rep(NA, 5),diff(tmin_mean, lag=5)),
                diff_tmin_6=c(rep(NA, 6), diff(tmin_mean, lag=6)),
                roll_tmin_2=RcppRoll::roll_meanr(tmin_mean, n=2, fill=NA, na.rm=T),
                roll_tmin_3=RcppRoll::roll_meanr(tmin_mean, n=3, fill=NA, na.rm=T),
                roll_tmin_4=RcppRoll::roll_meanr(tmin_mean, n=4, fill=NA, na.rm=T),
                roll_tmin_5=RcppRoll::roll_meanr(tmin_mean, n=5, fill=NA, na.rm=T),
                roll_tmin_6=RcppRoll::roll_meanr(tmin_mean, n=6, fill=NA, na.rm=T),
                lagg_prcp_1=lag(prcp_mean),
                diff_prcp_1=c(NA, diff(prcp_mean, lag=1)),
                diff_prcp_2=c(rep(NA, 2), diff(prcp_mean, lag=2)),
                diff_prcp_3=c(rep(NA, 3),diff(prcp_mean, lag=3)),
                diff_prcp_4=c(rep(NA, 4), diff(prcp_mean, lag=4)),
                diff_prcp_5=c(rep(NA, 5),diff(prcp_mean, lag=5)),
                diff_prcp_6=c(rep(NA, 6), diff(prcp_mean, lag=6)),
                roll_prcp_2=RcppRoll::roll_meanr(prcp_mean, n=2, fill=NA, na.rm=T),
                roll_prcp_3=RcppRoll::roll_meanr(prcp_mean, n=3, fill=NA, na.rm=T),
                roll_prcp_4=RcppRoll::roll_meanr(prcp_mean, n=4, fill=NA, na.rm=T),
                roll_prcp_5=RcppRoll::roll_meanr(prcp_mean, n=5, fill=NA, na.rm=T),
                roll_prcp_6=RcppRoll::roll_meanr(prcp_mean, n=6, fill=NA, na.rm=T),
                lagg_vp_1=lag(vp_mean),
                diff_vp_1=c(NA, diff(vp_mean, lag=1)),
                diff_vp_2=c(rep(NA, 2), diff(vp_mean, lag=2)),
                diff_vp_3=c(rep(NA, 3),diff(vp_mean, lag=3)),
                diff_vp_4=c(rep(NA, 4), diff(vp_mean, lag=4)),
                diff_vp_5=c(rep(NA, 5),diff(vp_mean, lag=5)),
                diff_vp_6=c(rep(NA, 6), diff(vp_mean, lag=6)),
                roll_vp_2=RcppRoll::roll_meanr(vp_mean, n=2, fill=NA, na.rm=T),
                roll_vp_3=RcppRoll::roll_meanr(vp_mean, n=3, fill=NA, na.rm=T),
                roll_vp_4=RcppRoll::roll_meanr(vp_mean, n=4, fill=NA, na.rm=T),
                roll_vp_5=RcppRoll::roll_meanr(vp_mean, n=5, fill=NA, na.rm=T),
                roll_vp_6=RcppRoll::roll_meanr(vp_mean, n=6, fill=NA, na.rm=T),
                lagg_srad_1=lag(srad_mean),
                diff_srad_1=c(NA, diff(srad_mean, lag=1)),
                diff_srad_2=c(rep(NA, 2), diff(srad_mean, lag=2)),
                diff_srad_3=c(rep(NA, 3),diff(srad_mean, lag=3)),
                diff_srad_4=c(rep(NA, 4), diff(srad_mean, lag=4)),
                diff_srad_5=c(rep(NA, 5),diff(srad_mean, lag=5)),
                diff_srad_6=c(rep(NA, 6), diff(srad_mean, lag=6)),
                roll_srad_2=RcppRoll::roll_meanr(srad_mean, n=2, fill=NA, na.rm=T),
                roll_srad_3=RcppRoll::roll_meanr(srad_mean, n=3, fill=NA, na.rm=T),
                roll_srad_4=RcppRoll::roll_meanr(srad_mean, n=4, fill=NA, na.rm=T),
                roll_srad_5=RcppRoll::roll_meanr(srad_mean, n=5, fill=NA, na.rm=T),
                roll_srad_6=RcppRoll::roll_meanr(srad_mean, n=6, fill=NA, na.rm=T),
                lagg_par_1=lag(par_tot),
                diff_par_1=c(NA, diff(par_tot, lag=1)),
                diff_par_2=c(rep(NA, 2), diff(par_tot, lag=2)),
                diff_par_3=c(rep(NA, 3),diff(par_tot, lag=3)),
                diff_par_4=c(rep(NA, 4), diff(par_tot, lag=4)),
                diff_par_5=c(rep(NA, 5),diff(par_tot, lag=5)),
                diff_par_6=c(rep(NA, 6), diff(par_tot, lag=6)),
                roll_par_2=RcppRoll::roll_sumr(par_tot, n=2, fill=NA, na.rm=T),
                roll_par_3=RcppRoll::roll_sumr(par_tot, n=3, fill=NA, na.rm=T),
                roll_par_4=RcppRoll::roll_sumr(par_tot, n=4, fill=NA, na.rm=T),
                roll_par_5=RcppRoll::roll_sumr(par_tot, n=5, fill=NA, na.rm=T),
                roll_par_6=RcppRoll::roll_sumr(par_tot, n=6, fill=NA, na.rm=T)
  )

daymet$key<-paste(daymet$date, daymet$pix)
daymet<-ungroup(daymet)
daymet<-dplyr::select(daymet, -date,  -pix)

write.csv(daymet, "output/daymet_monthly_climate_seagrant_plots.csv", row.names=F)

##All GCE data
files<-list.files("data/")
files<-files[grep("daymet_daily_climate_allgce", files)]
daymet<-read_csv(paste0("data/", files[1]))
daymet$date<-as.Date(daymet$date)
daymet$mo<-format(daymet$date, "%m")
daymet$yr<-format(daymet$date, "%Y")
daymet$pix<-paste(daymet$utm_east, daymet$utm_north)
daymet<-dplyr::select(daymet, -'system:index', -utm_east, -utm_north, -elevation, -.geo, -date, -swe)
daymet<-daymet %>% group_by(yr, mo, pix) %>% summarise_all(funs(mean, sum))
daymet$date<-as.Date(paste(daymet$yr, daymet$mo, "15", sep="-"))
daymet<-ungroup(daymet)
daymet<-dplyr::select(daymet, -yr, -mo)
daymet$par_sum<-2.114*daymet$srad_sum
daymet$par_tot<-daymet$dayl_sum*daymet$par_sum/1000000

daymet<-arrange(daymet, pix, date)
daymet<-daymet %>% dplyr::group_by( pix) %>%
  dplyr::mutate(lagg_tmax_1=lag(tmax_mean),
                diff_tmax_1=c(NA, diff(tmax_mean, lag=1)),
                diff_tmax_2=c(rep(NA, 2), diff(tmax_mean, lag=2)),
                diff_tmax_3=c(rep(NA, 3),diff(tmax_mean, lag=3)),
                diff_tmax_4=c(rep(NA, 4), diff(tmax_mean, lag=4)),
                diff_tmax_5=c(rep(NA, 5),diff(tmax_mean, lag=5)),
                diff_tmax_6=c(rep(NA, 6), diff(tmax_mean, lag=6)),
                roll_tmax_2=RcppRoll::roll_meanr(tmax_mean, n=2, fill=NA, na.rm=T),
                roll_tmax_3=RcppRoll::roll_meanr(tmax_mean, n=3, fill=NA, na.rm=T),
                roll_tmax_4=RcppRoll::roll_meanr(tmax_mean, n=4, fill=NA, na.rm=T),
                roll_tmax_5=RcppRoll::roll_meanr(tmax_mean, n=5, fill=NA, na.rm=T),
                roll_tmax_6=RcppRoll::roll_meanr(tmax_mean, n=6, fill=NA, na.rm=T),
                lagg_tmin_1=lag(tmin_mean),
                diff_tmin_1=c(NA, diff(tmin_mean, lag=1)),
                diff_tmin_2=c(rep(NA, 2), diff(tmin_mean, lag=2)),
                diff_tmin_3=c(rep(NA, 3),diff(tmin_mean, lag=3)),
                diff_tmin_4=c(rep(NA, 4), diff(tmin_mean, lag=4)),
                diff_tmin_5=c(rep(NA, 5),diff(tmin_mean, lag=5)),
                diff_tmin_6=c(rep(NA, 6), diff(tmin_mean, lag=6)),
                roll_tmin_2=RcppRoll::roll_meanr(tmin_mean, n=2, fill=NA, na.rm=T),
                roll_tmin_3=RcppRoll::roll_meanr(tmin_mean, n=3, fill=NA, na.rm=T),
                roll_tmin_4=RcppRoll::roll_meanr(tmin_mean, n=4, fill=NA, na.rm=T),
                roll_tmin_5=RcppRoll::roll_meanr(tmin_mean, n=5, fill=NA, na.rm=T),
                roll_tmin_6=RcppRoll::roll_meanr(tmin_mean, n=6, fill=NA, na.rm=T),
                lagg_prcp_1=lag(prcp_mean),
                diff_prcp_1=c(NA, diff(prcp_mean, lag=1)),
                diff_prcp_2=c(rep(NA, 2), diff(prcp_mean, lag=2)),
                diff_prcp_3=c(rep(NA, 3),diff(prcp_mean, lag=3)),
                diff_prcp_4=c(rep(NA, 4), diff(prcp_mean, lag=4)),
                diff_prcp_5=c(rep(NA, 5),diff(prcp_mean, lag=5)),
                diff_prcp_6=c(rep(NA, 6), diff(prcp_mean, lag=6)),
                roll_prcp_2=RcppRoll::roll_meanr(prcp_mean, n=2, fill=NA, na.rm=T),
                roll_prcp_3=RcppRoll::roll_meanr(prcp_mean, n=3, fill=NA, na.rm=T),
                roll_prcp_4=RcppRoll::roll_meanr(prcp_mean, n=4, fill=NA, na.rm=T),
                roll_prcp_5=RcppRoll::roll_meanr(prcp_mean, n=5, fill=NA, na.rm=T),
                roll_prcp_6=RcppRoll::roll_meanr(prcp_mean, n=6, fill=NA, na.rm=T),
                lagg_vp_1=lag(vp_mean),
                diff_vp_1=c(NA, diff(vp_mean, lag=1)),
                diff_vp_2=c(rep(NA, 2), diff(vp_mean, lag=2)),
                diff_vp_3=c(rep(NA, 3),diff(vp_mean, lag=3)),
                diff_vp_4=c(rep(NA, 4), diff(vp_mean, lag=4)),
                diff_vp_5=c(rep(NA, 5),diff(vp_mean, lag=5)),
                diff_vp_6=c(rep(NA, 6), diff(vp_mean, lag=6)),
                roll_vp_2=RcppRoll::roll_meanr(vp_mean, n=2, fill=NA, na.rm=T),
                roll_vp_3=RcppRoll::roll_meanr(vp_mean, n=3, fill=NA, na.rm=T),
                roll_vp_4=RcppRoll::roll_meanr(vp_mean, n=4, fill=NA, na.rm=T),
                roll_vp_5=RcppRoll::roll_meanr(vp_mean, n=5, fill=NA, na.rm=T),
                roll_vp_6=RcppRoll::roll_meanr(vp_mean, n=6, fill=NA, na.rm=T),
                lagg_srad_1=lag(srad_mean),
                diff_srad_1=c(NA, diff(srad_mean, lag=1)),
                diff_srad_2=c(rep(NA, 2), diff(srad_mean, lag=2)),
                diff_srad_3=c(rep(NA, 3),diff(srad_mean, lag=3)),
                diff_srad_4=c(rep(NA, 4), diff(srad_mean, lag=4)),
                diff_srad_5=c(rep(NA, 5),diff(srad_mean, lag=5)),
                diff_srad_6=c(rep(NA, 6), diff(srad_mean, lag=6)),
                roll_srad_2=RcppRoll::roll_meanr(srad_mean, n=2, fill=NA, na.rm=T),
                roll_srad_3=RcppRoll::roll_meanr(srad_mean, n=3, fill=NA, na.rm=T),
                roll_srad_4=RcppRoll::roll_meanr(srad_mean, n=4, fill=NA, na.rm=T),
                roll_srad_5=RcppRoll::roll_meanr(srad_mean, n=5, fill=NA, na.rm=T),
                roll_srad_6=RcppRoll::roll_meanr(srad_mean, n=6, fill=NA, na.rm=T),
                lagg_par_1=lag(par_tot),
                diff_par_1=c(NA, diff(par_tot, lag=1)),
                diff_par_2=c(rep(NA, 2), diff(par_tot, lag=2)),
                diff_par_3=c(rep(NA, 3),diff(par_tot, lag=3)),
                diff_par_4=c(rep(NA, 4), diff(par_tot, lag=4)),
                diff_par_5=c(rep(NA, 5),diff(par_tot, lag=5)),
                diff_par_6=c(rep(NA, 6), diff(par_tot, lag=6)),
                roll_par_2=RcppRoll::roll_sumr(par_tot, n=2, fill=NA, na.rm=T),
                roll_par_3=RcppRoll::roll_sumr(par_tot, n=3, fill=NA, na.rm=T),
                roll_par_4=RcppRoll::roll_sumr(par_tot, n=4, fill=NA, na.rm=T),
                roll_par_5=RcppRoll::roll_sumr(par_tot, n=5, fill=NA, na.rm=T),
                roll_par_6=RcppRoll::roll_sumr(par_tot, n=6, fill=NA, na.rm=T)
  )

daymet$key<-paste(daymet$date, daymet$pix)
daymet<-ungroup(daymet)
daymet<-dplyr::select(daymet, -date,  -pix)

#i<-2
## interate through the rest
for( i in 2:length(files)){
  dat<-read_csv(paste0("data/", files[i]))
  dat$date<-as.Date(dat$date)
  dat$mo<-format(dat$date, "%m")
  dat$yr<-format(dat$date, "%Y")
  dat$pix<-paste(dat$utm_east, dat$utm_north)
  dat<-dplyr::select(dat, -'system:index', -utm_east, -utm_north, -elevation, -.geo, -date, -swe)
  dat<-dat %>% group_by(yr, mo, pix) %>% summarise_all(funs(mean, sum))
  dat$date<-as.Date(paste(dat$yr, dat$mo, "15", sep="-"))
  dat<-ungroup(dat)
  dat<-dplyr::select(dat, -yr, -mo)
  dat$par_sum<-2.114*dat$srad_sum
  dat$par_tot<-dat$dayl_sum*dat$par_sum/1000000
  
  dat<-arrange(dat, pix, date)
  dat<-dat %>% dplyr::group_by( pix) %>%
    dplyr::mutate(lagg_tmax_1=lag(tmax_mean),
                  diff_tmax_1=c(NA, diff(tmax_mean, lag=1)),
                  diff_tmax_2=c(rep(NA, 2), diff(tmax_mean, lag=2)),
                  diff_tmax_3=c(rep(NA, 3),diff(tmax_mean, lag=3)),
                  diff_tmax_4=c(rep(NA, 4), diff(tmax_mean, lag=4)),
                  diff_tmax_5=c(rep(NA, 5),diff(tmax_mean, lag=5)),
                  diff_tmax_6=c(rep(NA, 6), diff(tmax_mean, lag=6)),
                  roll_tmax_2=RcppRoll::roll_meanr(tmax_mean, n=2, fill=NA, na.rm=T),
                  roll_tmax_3=RcppRoll::roll_meanr(tmax_mean, n=3, fill=NA, na.rm=T),
                  roll_tmax_4=RcppRoll::roll_meanr(tmax_mean, n=4, fill=NA, na.rm=T),
                  roll_tmax_5=RcppRoll::roll_meanr(tmax_mean, n=5, fill=NA, na.rm=T),
                  roll_tmax_6=RcppRoll::roll_meanr(tmax_mean, n=6, fill=NA, na.rm=T),
                  lagg_tmin_1=lag(tmin_mean),
                  diff_tmin_1=c(NA, diff(tmin_mean, lag=1)),
                  diff_tmin_2=c(rep(NA, 2), diff(tmin_mean, lag=2)),
                  diff_tmin_3=c(rep(NA, 3),diff(tmin_mean, lag=3)),
                  diff_tmin_4=c(rep(NA, 4), diff(tmin_mean, lag=4)),
                  diff_tmin_5=c(rep(NA, 5),diff(tmin_mean, lag=5)),
                  diff_tmin_6=c(rep(NA, 6), diff(tmin_mean, lag=6)),
                  roll_tmin_2=RcppRoll::roll_meanr(tmin_mean, n=2, fill=NA, na.rm=T),
                  roll_tmin_3=RcppRoll::roll_meanr(tmin_mean, n=3, fill=NA, na.rm=T),
                  roll_tmin_4=RcppRoll::roll_meanr(tmin_mean, n=4, fill=NA, na.rm=T),
                  roll_tmin_5=RcppRoll::roll_meanr(tmin_mean, n=5, fill=NA, na.rm=T),
                  roll_tmin_6=RcppRoll::roll_meanr(tmin_mean, n=6, fill=NA, na.rm=T),
                  lagg_prcp_1=lag(prcp_mean),
                  diff_prcp_1=c(NA, diff(prcp_mean, lag=1)),
                  diff_prcp_2=c(rep(NA, 2), diff(prcp_mean, lag=2)),
                  diff_prcp_3=c(rep(NA, 3),diff(prcp_mean, lag=3)),
                  diff_prcp_4=c(rep(NA, 4), diff(prcp_mean, lag=4)),
                  diff_prcp_5=c(rep(NA, 5),diff(prcp_mean, lag=5)),
                  diff_prcp_6=c(rep(NA, 6), diff(prcp_mean, lag=6)),
                  roll_prcp_2=RcppRoll::roll_meanr(prcp_mean, n=2, fill=NA, na.rm=T),
                  roll_prcp_3=RcppRoll::roll_meanr(prcp_mean, n=3, fill=NA, na.rm=T),
                  roll_prcp_4=RcppRoll::roll_meanr(prcp_mean, n=4, fill=NA, na.rm=T),
                  roll_prcp_5=RcppRoll::roll_meanr(prcp_mean, n=5, fill=NA, na.rm=T),
                  roll_prcp_6=RcppRoll::roll_meanr(prcp_mean, n=6, fill=NA, na.rm=T),
                  lagg_vp_1=lag(vp_mean),
                  diff_vp_1=c(NA, diff(vp_mean, lag=1)),
                  diff_vp_2=c(rep(NA, 2), diff(vp_mean, lag=2)),
                  diff_vp_3=c(rep(NA, 3),diff(vp_mean, lag=3)),
                  diff_vp_4=c(rep(NA, 4), diff(vp_mean, lag=4)),
                  diff_vp_5=c(rep(NA, 5),diff(vp_mean, lag=5)),
                  diff_vp_6=c(rep(NA, 6), diff(vp_mean, lag=6)),
                  roll_vp_2=RcppRoll::roll_meanr(vp_mean, n=2, fill=NA, na.rm=T),
                  roll_vp_3=RcppRoll::roll_meanr(vp_mean, n=3, fill=NA, na.rm=T),
                  roll_vp_4=RcppRoll::roll_meanr(vp_mean, n=4, fill=NA, na.rm=T),
                  roll_vp_5=RcppRoll::roll_meanr(vp_mean, n=5, fill=NA, na.rm=T),
                  roll_vp_6=RcppRoll::roll_meanr(vp_mean, n=6, fill=NA, na.rm=T),
                  lagg_srad_1=lag(srad_mean),
                  diff_srad_1=c(NA, diff(srad_mean, lag=1)),
                  diff_srad_2=c(rep(NA, 2), diff(srad_mean, lag=2)),
                  diff_srad_3=c(rep(NA, 3),diff(srad_mean, lag=3)),
                  diff_srad_4=c(rep(NA, 4), diff(srad_mean, lag=4)),
                  diff_srad_5=c(rep(NA, 5),diff(srad_mean, lag=5)),
                  diff_srad_6=c(rep(NA, 6), diff(srad_mean, lag=6)),
                  roll_srad_2=RcppRoll::roll_meanr(srad_mean, n=2, fill=NA, na.rm=T),
                  roll_srad_3=RcppRoll::roll_meanr(srad_mean, n=3, fill=NA, na.rm=T),
                  roll_srad_4=RcppRoll::roll_meanr(srad_mean, n=4, fill=NA, na.rm=T),
                  roll_srad_5=RcppRoll::roll_meanr(srad_mean, n=5, fill=NA, na.rm=T),
                  roll_srad_6=RcppRoll::roll_meanr(srad_mean, n=6, fill=NA, na.rm=T), 
                  lagg_par_1=lag(par_tot),
                  diff_par_1=c(NA, diff(par_tot, lag=1)),
                  diff_par_2=c(rep(NA, 2), diff(par_tot, lag=2)),
                  diff_par_3=c(rep(NA, 3),diff(par_tot, lag=3)),
                  diff_par_4=c(rep(NA, 4), diff(par_tot, lag=4)),
                  diff_par_5=c(rep(NA, 5),diff(par_tot, lag=5)),
                  diff_par_6=c(rep(NA, 6), diff(par_tot, lag=6)),
                  roll_par_2=RcppRoll::roll_sumr(par_tot, n=2, fill=NA, na.rm=T),
                  roll_par_3=RcppRoll::roll_sumr(par_tot, n=3, fill=NA, na.rm=T),
                  roll_par_4=RcppRoll::roll_sumr(par_tot, n=4, fill=NA, na.rm=T),
                  roll_par_5=RcppRoll::roll_sumr(par_tot, n=5, fill=NA, na.rm=T),
                  roll_par_6=RcppRoll::roll_sumr(par_tot, n=6, fill=NA, na.rm=T)
    )
  
  dat$key<-paste(dat$date, dat$pix)
  dat<-ungroup(dat)
  dat<-dplyr::select(dat, -date,  -pix)
  daymet<-rbind(daymet, dat)
}

write.csv(daymet,"output/daymet_monthly_climate_newpixels.csv", row.names=F)
daymet$date<-str_split(daymet$key, " ", simplify = T)[,1]
daymet$year<-format(as.Date(daymet$date), "%Y")
table(daymet$year)