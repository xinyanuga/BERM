###This is the cluster version!!!
###A function to calculate soil temperature from Land surface temperature and 
##air temperature and estimate total soil growing degrees and the date fo spring green-up
library(tidyverse); library(mgcv); library(RcppRoll)

setup_heat<-function(data, gce=F){
  heat.l8<-data
  if(gce==F){
    heat.l8<-ungroup(heat.l8)
    heat.l8<-dplyr::select(heat.l8,  date, mo, year, utm_east, utm_north, pix, lst, elevation)
    heat.l8$mo<-as.numeric(heat.l8$mo)
    heat.l8<-heat.l8 %>% dplyr::group_by(date, mo, year, pix) %>% 
      dplyr::summarise( utm_east=mean(utm_east, na.rm=T), utm_north=mean(utm_north, na.rm=T),lst=mean(lst, na.rm=T), elevation=mean(elevation))
  }
  
  soil.mod<-readRDS("output/soiltemp_from_lst_and_air_lags.rds")
  
  weeks<-read.csv("data/week_lookup_table.csv")
  climate<-read.csv("data/nws_sapelo_daily_jan2013_sept2020.csv")
  climate$DATE<-as.Date(climate$DATE)
  climate<-dplyr::select(climate, DATE,DAPR,MDPR,PRCP,TMAX,TMIN)
  names(climate)<-c("date", "dapr", "mdpr", "prcp", "tmax", "tmin")
  climate$mo<-format(climate$date, "%m")
  
  ## a function to identify short runs of missing values in a time series
  seqle <- function(x,incr=1) { 
    if(!is.numeric(x)) x <- as.numeric(x) 
    n <- length(x)  
    y <- x[-1L] != x[-n] + incr 
    i <- c(which(y|is.na(y)),n) 
    list(lengths = diff(c(0L,i)),
         values = x[head(c(0L,i)+1L,-1L)]) 
  } 
  
  ##fill in short runs of missing values in the air temperature data with surrounding means
  ##see which values are NA, apply seqle to see if they are in sequence, then fill them in if they're short runs
  x<-which(is.na(climate$tmin))
  s<-seqle(x)
  mean<-roll_mean(climate$tmin, n = 9, by = 1,  align = "center", fill=NA, na.rm=T)
  for(i in seq_along(s$lengths)){
    if(s$lengths[i]<4){
      n<-seq(s$values[i], (s$values[i]+s$lengths[i]-1))
      for(j in seq_along(n)){
        climate$tmin[n]<-mean[n] }
    }
  }
  x<-which(is.na(climate$tmin))
  s<-seqle(x)
  
  x<-which(is.na(climate$tmax))
  s<-seqle(x)
  mean<-roll_mean(climate$tmax, n = 9, by = 1,  align = "center", fill=NA, na.rm=T)
  for(i in seq_along(s$lengths)){
    if(s$lengths[i]<4){
      n<-seq(s$values[i], (s$values[i]+s$lengths[i]-1))
      for(j in seq_along(n)){
        climate$tmax[n]<-mean[n] }
    }
  }
  x<-which(is.na(climate$tmax))
  s<-seqle(x)
  
  ##calculate lags in air temperature, so that temps from past days can influence today soil temp
  climate<-arrange(climate, date)
  climate$tmean<-(climate$tmax+climate$tmin)/2
  climate$ran<-climate$tmax-climate$tmin
  climate<-dplyr::select(climate, date, tmean, ran)
  climate$lag1=lag(as.numeric(climate$tmean))
  climate$lag2<-lag(climate$lag1)
  climate$lag3<-lag(climate$lag2)
  climate$lag4<-lag(climate$lag3)
  climate$lag5<-lag(climate$lag4)
  climate$lag6<-lag(climate$lag5)
  climate$lag7<-lag(climate$lag6)
  climate$lag8<-lag(climate$lag7)
  climate$lag9<-lag(climate$lag8)
  climate$lag10<-lag(climate$lag9)
  climate$lag11<-lag(climate$lag10)
  elevation<-unique(heat.l8[,c("pix", "elevation", "utm_east", "utm_north")])
  
  ##merge the air temp climate data and the pixel data, along with the look up table that has the average differences between air temp and LST by pixel to gap fill missing values
  diffs<-merge(climate, heat.l8, by="date")
  weeks$date<-paste(weeks$year, weeks$month, weeks$day, sep="-")
  diffs$key<-findInterval(as.Date(diffs$date), as.Date(weeks$date))
  wks<-dplyr::select(weeks, biweek, key)
  diffs<-merge(diffs, wks, by="key")
  diffs$diff<-diffs$tmean-diffs$lst
  diffs<-diffs %>% dplyr::group_by(pix, biweek) %>% dplyr::summarise(diff=mean(diff, na.rm=T))
  diffs<-ungroup(diffs)
  diffs$key<-paste(diffs$pix, diffs$biweek)
  diffs<-dplyr::select(diffs, -pix, -biweek)
  
  heat.l8<-heat.l8[!is.na(heat.l8$date),]
  pix<-unique(heat.l8$pix)
  wks<-data.frame(pix=rep(pix, each=nrow(weeks)), week_date=rep(weeks$date,length(pix)),
                  biweek=rep(weeks$biweek, length(pix)), key=rep(weeks$key, length(pix)))
  wks$key<-paste(wks$pix, wks$key)
  wks<-arrange(wks, pix,week_date)
  wks<-dplyr::select(wks, -pix)
  heat.l8$key<-findInterval(as.Date(heat.l8$date), as.Date(weeks$date))
  heat.l8$key<-paste(heat.l8$pix, heat.l8$key)
  heat.l8<-merge(heat.l8, wks, by="key", all.x=T, all.y=T)
  heat.l8$date<-ifelse(is.na(heat.l8$date), as.Date(heat.l8$week_date)-7, heat.l8$date)
  heat.l8$date<-as.Date(heat.l8$date, origin="1970-01-01")
  heat.l8$mo<-format(heat.l8$date, "%m")
  heat.l8$year<-format(heat.l8$date, "%Y")
  heat.l8<-arrange(heat.l8, pix,year, biweek)
  if(gce==T){
    x<-str_split(heat.l8$key, " ", simplify = T)
    heat.l8$utm_east<- x[,1]
    heat.l8$utm_north<-x[,2]
    heat.l8$pix<- paste(heat.l8$utm_east, heat.l8$utm_north)
  }
  if(gce==F){
    x<-str_split(heat.l8$key, " ", simplify = T)
    heat.l8$pix<- x[,1]
  }
  heat.l8$key<-paste(heat.l8$pix, heat.l8$biweek)
  heat.l8<-merge(heat.l8, climate, by="date")
  heat.l8<-full_join(heat.l8,diffs, by="key",all.y=TRUE, all.x=T)
  heat.l8$doy<-as.numeric(format(heat.l8$date, "%j"))
  heat.l8<-dplyr::select(heat.l8, -elevation, -utm_east, -utm_north)
  heat.l8<-left_join(heat.l8, elevation, by="pix")
  
  ##gap fill LST as the difference beteen mean air temp and the average weekly difference between air temp and LST for that pixel
  ##also gap fill unreasonable values (diff between air and lst >12 deg Celsius)
  heat.l8$lst<-ifelse(is.na(heat.l8$lst), heat.l8$tmean-heat.l8$diff, heat.l8$lst)
  heat.l8$lst<-ifelse(abs(heat.l8$lst-heat.l8$tmean)>12, heat.l8$tmean-heat.l8$diff, heat.l8$lst)
  
  
  ##Now that LST is gap filled, interpolate to a daily estimate
  ##dates to predict over
  dates<-seq.Date(as.Date("2013-03-28"), as.Date(max(heat.l8$date)), by= "1 day")
  
  ##predict daily estimate through approx, but much faster than for loop
  ##this allows dplyr to return a result that is greater > 1 per group
  ##the approx implementation just makes a line between adjacent obs, no local weighting
  heat.l8<-heat.l8[!is.na(heat.l8$lst),]
  heat.l8<-dplyr::arrange(heat.l8,  date, pix)
  heat.l8$pixelev<-heat.l8$elevation
  out<-heat.l8 %>% dplyr::group_by(pix, pixelev) %>% dplyr::filter(length(date)>3) %>%
    do(data.frame(date=dates,
                  lst= approx(as.numeric(.$date),  .$lst, 
                              xout=as.numeric(dates), method="linear")$y
    )
    )
  
  
  out$date<-as.Date(out$date, origin="1970-01-01")
  out$year<-as.numeric(format(out$date, "%Y"))
  out$mo<-as.numeric(format(out$date, "%m"))
  out$doy<-as.numeric(format(out$date, "%j"))
  
  out<-merge(out, climate, by="date")
  out$pixelev.x<-out$pixelev
  
  ##estimate the soil temperature from the soil temperature model previously created, which depends on LST and lagged air temp
  out$soil.temp<-(predict(soil.mod, newdata=out))
  
  ##calculate the annual date of greenup on a pixelwise basis
  ##the base temperature identified in O'Connell et al, 2021 Ecosystems
  base<-9.9
  ##calculate total soil growing degree days
  out$dd.soil<-ifelse(out$soil.temp<base,0,out$soil.temp-base)
  out$dd.lst<-ifelse(out$lst<base,0,out$lst-base)
  out$mo<-format(out$date, "%m")
  
  ##estimate the date of green up as the date that the sums of total growing degree days exceeds 202;
  years<-unique(out$year); years<-years[years>2013]
  pix<-as.character(unique(out$pix))
  tdd<-data.frame(pix=rep("A", length(years)*length(pix)), year=rep(9999, length(years)*length(pix)),
                  greenup.lst=rep(as.Date("1970-01-01"), length(years)*length(pix)),
                  greenup.soil=rep(as.Date("1970-01-01"), length(years)*length(pix)), stringsAsFactors = F)
  n<-1
  for (i in 1:length(years)){
    for (j in 1:length(pix)){
      start<-as.Date(paste0(years[i], "-01-01"), origin="1970-01-01")
      end<-as.Date(paste0(years[i], "-05-15"), origin="1970-01-01")
      temp<-out[out$date>=start&out$date<=end&out$pix==pix[j],]
      tdd.poss.soil<-cumsum(temp$dd.soil)
      tdd.poss.lst<-cumsum(temp$dd.lst)
      
      tdd$pix[n]<-pix[j]
      tdd$year[n]<-years[i]
      tdd$greenup.soil[n]<-temp$date[min(which(tdd.poss.soil>202))]
      tdd$greenup.lst[n]<-temp$date[min(which(tdd.poss.lst>202))]
      n<-1+n
    }
  }
  
  ##create averages of observed winter lst, in case we need it for later model estimates
  heat.l8$year<-as.numeric(heat.l8$year)
  heat2.l8<- out %>% dplyr::group_by(pix,year) %>% 
    dplyr::filter (mo %in% c("12","01","02")) 
  heat2.l8$year<-ifelse(heat2.l8$mo=="12", heat2.l8$year+1, heat2.l8$year)
  heat.l8<- heat2.l8 %>% dplyr::group_by(pix, year) %>% 
    dplyr::summarise(winlst.sum=sum(lst, na.rm=T),
                     winlst.avg=mean(lst, na.rm=T))
  
  heat2.l8<- out %>% dplyr::group_by(pix,year) %>% 
    dplyr::filter(as.numeric(format(date, "%j"))<50)  %>% 
    dplyr::summarise(lst2.sum=sum(lst),lst2.avg=mean(lst, na.rm=T))
  heat.l8<-full_join(heat.l8, heat2.l8)
  heat.l8<-ungroup(heat.l8)
  heat.l8<-left_join(heat.l8, tdd)
  
  return(heat.l8)
  
}
