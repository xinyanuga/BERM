##need "flood" model from landsat is flooring.r in memory before running this script
library(tidyverse)

##Set the working directory; change path to your working directory
setwd("/home/jessica/UGA/Projects/BERM")

#Load the field data
veg<-read.csv("data/seagrant_core_vegplot_combined.csv", header=T)
#veg<-merge(veg, plots, by="plot")
veg<-arrange(veg, date)
veg$date<-as.Date(veg$date)
veg$year<-as.numeric(format(veg$date, "%Y"))
veg$mo<-as.numeric(format(veg$date, "%m"))
veg$junc_live_cover<-ifelse(veg$site=="fluxa",0,veg$junc_live_cover)
veg$brown<-ifelse(veg$brown=="yes", 1, ifelse(veg$brown=="iron precip, green", 0.5, 0))
veg$brown[which(is.na(veg$brown))]<-0

##smooth veg data series to estimate info on all valid landsat dates
#i<-1; j<-1; k<-1
var<-c("percentN", "bgm2.areascaled","agm2.core.areascaled", "bgm2.core.stemscaled",
       "agm2.core.stemscaled", "agm2.allom", "bgm2.allometric.rootshoot",
       "bgm2.allometric.rootgreenshoot", "agm2.allom.new", "bgm2.allo.new.rootshoot",
       "bgm2.allo.new.rootgreenshoot",  "rhizm2.allometric.rootshoot",
       "rhizm2.allometric.rootgreenshoot", "rhizm2.allo.new.rootshoot",
       "rhizm2.allo.new.rootgreenshoot","rhizm2.areascaled","rhizm2.core.stemscaled",
       "chlorophyll", "LAI")
loc<-c("fb1", "fb2", "fb3", "fb4", "fb5", "fb6", "fb7", "fb8", "fb9",
       "m1", "m2", "m3", "m4", "m5", "m6", "s1", "s2", "s3", "s4", "s5", "s6",
       "sk1", "sk2", "sk3", "sk4", "sk5", "sk6", "sk7", "sk8", "sk9",
       "t1",  "t2", "t3", "t4", "t5", "t6","t7", "t8",
       "ug1", "ug2", "ug3", "ug4", "ug5", "ug6","ug7", "ug8", "ug9")
dates<-seq.Date(as.Date("2013-05-01"), as.Date("2020-01-01"), by= "1 week")
#dates<-dates[as.numeric(format(dates, "%m"))>4&as.numeric(format(dates, "%m"))<11]
out<-data.frame(matrix(ncol=length(var)+2, nrow=5000))
names(out)<-c("date", "plot", var)
loop<-1

for (j in seq_along(loc)){
  for (i in seq_along(var)){
    x<-veg[veg$plot==loc[j], c("date", var[i])]
    x$date<-as.numeric(x$date)
    x<-x[complete.cases(x),]
    
    if(nrow(x)>1){
      smoothed <- approx(x$date,x[,var[i]], xout=as.numeric(dates),
                         method="linear"
      )
      smoothy<-smoothed}
    num<-length(smoothy$x)
    out[loop:(loop+num-1),"plot"]<-rep(loc[j], num)
    out[loop:(loop+num-1),"date"]<-as.Date(smoothy$x, origin = "1970-01-01")
    out[loop:(loop+num-1), var[i]]<-(smoothy$y)
    
  }
  loop<-loop+num
}
out$date<-as.Date(out$date, origin="1970-01-01")
out$year<-as.numeric(format(out$date, "%Y"))
out$mo<-as.numeric(format(out$date, "%m"))
out$site<-substr(out$plot, 1,2)
out$site[!(out$site %in% c("fb", "sk", "ug"))]<- "fa"

##remove unsampled years at other sites
out<-out[out$site=="fa"|(out$site=="fb"&out$year==2016)|(out$site=="ug"&out$year==2016)|(out$site=="sk"&out$year==2016),]
out[out$plot=="t8"&out$year|out$plot=="t7","percentN"]<-NA
out[out$plot=="t8"&out$year|out$plot=="t7","chlorophyll"]<-NA
out[out$plot=="t8"&out$year|out$plot=="t7","LAI"]<-NA

out$chlorophyll[out$year<2016]<-NA
out$chlorophyll[out$year==2016&out$mo<5]<-NA
out$LAI[out$year==2016&out$mo<5]<-NA
out$LAI[out$year<2016]<-NA
out<-out[as.Date(out$date)>"2013-06-11",]
out<-out[!is.na(out$date),]

var<-c("rhizm2.allometric.rootshoot",
       "rhizm2.allometric.rootgreenshoot", "rhizm2.allo.new.rootshoot",
       "rhizm2.allo.new.rootgreenshoot","rhizm2.areascaled","rhizm2.core.stemscaled"
)

out[out$date<"2016-05-03",var]<-NA
out[out$date>"2016-10-29"&out$date<"2017-06-30",var]<-NA
out[out$date>"2017-07-18"&out$date<"2017-10-03",var]<-NA

var<-"percentN"
out[out$date<"2014-04-24",var]<-NA
out[out$date>"2014-10-03"&out$date<"2015-04-14",var]<-NA
out[out$date>"2015-09-23"&out$date<"2016-05-03",var]<-NA
out[out$date>"2015-09-23"&out$date<"2016-05-03",var]<-NA
out<-out[!is.na(out$agm2.allom), ]
out$site<-ifelse(out$site=="fb", "fluxb", out$site)
out$site<-ifelse(out$site=="fa", "fluxa", out$site)
out$site<-ifelse(out$site=="ug", "ugami", out$site)
out$site<-ifelse(out$site=="sk", "skida", out$site)

write.csv(out, "output/smoothed_veg_timeseries_weekly.csv", row.names = F)