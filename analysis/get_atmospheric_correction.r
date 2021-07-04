###A script to get atmospheric correction parameters from Julia Barsi's NASA web tool
###These parameters help convert top of atmosphere brightness temperature (B10 in Landsat 8, and B6 in Landsat 5)
###into surface temperature; USGS has released a provisional surface temp product 
###that will help us skip this step in the future
###This script is slow to run and, important, !CANNOT! be speed up
###Julia Barsi has asked users not to retrieve info more quickly than once every 2 mins, 
###a wait time is built into the download function and handles this for you, 
###so feel free to set your usage parameters, run it and walk away
###
###USAGE:Input file should be an Earth Explorer Landsat !cloud-free! data file with the .geo (for lat/long) and date fields retained
###      Change base directory and input file name to your file, 
###      change append.out to True or False as needed: 
###           TRUE means your running this after your first try and have old results from a previous try that you want to append to
###           Feel free to combine training, validation, etc observations from each 
###           landsat mission into 1 atmospheric parameter file, but you should separate the missions into different files
###           TRUE writes over your old file with a new combined file, 
###           it's a good idea to back-up your old file
###      change mission to 8 or 5 for Landsat 8 or 5
###      Change email to your email address, but don't give this script to others with your email address saved, 
###           Dr. Barsi will contact you via email if you misuse this tool
###      All other parameters are calcaulted automatically from the Landsat Earth Engine input file
###      After setting everything in the USE INPUT REQUIRED section below, save and run the script
###      Assumes that there is a project folder to hold everything (the "base directory") with the subfolders:
###                   "data", where the earth engine landsat data is held, and 
###                   "output", where the results of this script will go


library(httr); library(tidyverse)

##add new results to an old file?
append.out<-T
basedir<-"/home/jessica/UGA/Projects/seagrant/"

###files to read in####
##landsat, csv pts from google earth engine
lsat<-read_csv(paste0(basedir, "data/L8_allgce_2013_2019_noclouds.csv"))

##old data we want to append to
if (append.out==T){
  old<-read_csv(paste0(basedir, "output/atmospheric_correction_varsL8_gce.csv"))
  }

##create function to get atmospheric correction parameters from web
atmos<-function(date,lat, long, email, time= "16:00", mission="8", profile="2"){
    ##download atmospheric variables as calculated by web tool: atmcorr.gsfc.nasa.gov (see Barsi et al. 2003; Barsi et al. 2005 )
    ##In a for loop of dates, this takes AWHILE to run, there also is a deliberate 2 min wait at the end of the function to prevent overloading the webtool with requests, be patient
    ##The variables, unless local data are entered, were calculated as an average of 1 deg lat and long (31-32 degrees lat for example) 
    ##results were very close to measured conditions at the sapelo flux tower for the days I tested
    ##profile_option = 1  Use atmospheric profile for closest integer lat/long, 2 = Use interpolated atmospheric profile for given lat/long, I used 2 (see previous sentence)
    ##stdatm_option value="1" Use mid-latitude summer standard atmosphere for upper atmospheric profile , 2 =Use mid-latitude winter standard atmosphere for upper atmospheric profile
    ##  this function uses the date to decide which stdatm to use, doy<91 &doy > 273 =winter
    ##L57_option = which landsat, put value= 8, 7,  5, or 11  where the number is the landsat mission and 11 means Output only atmospheric profile, do not calculate effective radiances
    ##time "16:00" is approx landsat passover time on the US East coast, which varies between 15:55 and 16:05 GMT, you can verify this in the landsat meta-data
    
    ##transform the variables into the format expected by the web form
    date<-as.Date(date)
    yr<-as.character(format(date, "%Y"))
    mo<-as.character(format(date, "%m"))
    day<-as.character(format(date, "%d"))
    hr<- substr(time, 1,2)
    min<- substr(time, 4,5)
    
    ##use date to decide which atomsphereic profile to use, winter or summer
    stdatm<-ifelse(as.numeric(format(date, "%j"))<91, "1","2")
    stdatm<-ifelse(as.numeric(format(date, "%j"))>273, "1",stdatm)
    
    ##check first that at least mission was entered correctly
    if(!(as.numeric(mission) %in% c(8,7,5,11))){
        print("mission must be 8,7,5 or 11")
    } else{
        ##gather form variables needed on the web form into a list
        formvars<-list(
            "year" = yr,
            "month" = mo,
            "day" = day,
            "hour" = as.character(hr),
            "minute" = as.character(min),
            "L57_option" = as.character(mission),
            "thelat" = as.character(round(lat,4)),
            "thelong" = as.character(round(long,4)),
            "profile_option" = as.character(profile),
            "user_email" = as.character(email),
            "stdatm_option" = "1", 
            "submit" ="Calculate"
        )
        ##use POST to submit the form and retrieve the result from the web
        res <- POST("https://atmcorr.gsfc.nasa.gov/cgi-bin/atm_corr.pl", body = formvars, verbose())
        x<-content(res, "text")
        
        ##process the result into a data frame of atomospheric variables that this function will return
        y<-str_split(x, "<br>")
        y<-y[[1]]
        Iu<-y[grep("upwelling radiance",y)];  Iu<-str_split(Iu, ":")[[1]][2];Iu<-as.numeric(str_split(Iu, " W")[[1]][1])
        Id<-y[grep("downwelling radiance",y)];Id<-str_split(Id, ":")[[1]][2];Id<-as.numeric(str_split(Id, " W")[[1]][1])
        trans<-y[grep("atmospheric transmission",y)]; trans<-str_split(trans, ":")[[1]][2];trans<-as.numeric(str_split(trans, " W")[[1]][1])
        out<-data.frame("date"=date, "Iu"=Iu, "Id"=Id, "trans"=trans)
        
        ##add in an automatic wait, requested by the creater of the web tool to prevent server overload
        Sys.sleep(2*60)
        return(out)
    }
}

##Process landsat data to get the dates and lat long to use, assumes from Google Earth Engine
lsat$time<-substr(lsat$date,12,19)
lsat$date<-as.Date(substr(lsat$date,1,10))
lsat$year<-as.numeric(format(lsat$date, "%Y"))
lsat$doy<-as.numeric(format(lsat$date, "%j"))
lsat$mo<-as.numeric(format(lsat$date, "%m"))

##no need to atmosphereically correct bad pixels, so first subset to good data
##see https://landsat.usgs.gov/landsat-surface-reflectance-quality-assessment
lsat$pixel_qa<-R.utils::intToBin(lsat$pixel_qa)
lsat$qa_good<-ifelse(str_sub(lsat$pixel_qa,-6,-1)=="000010",T,F) ##"000010" from the right means no fill, yes clear, no water, no cloud shadow
lsat$qa_good<-ifelse(str_sub(lsat$pixel_qa,-6,-1) %in% c("000010", "000100"),T,F) ##"000010" from the right means no fill, yes clear, no water, no cloud shadow, "000100" is water
lsat$wet<-ifelse(str_sub(lsat$pixel_qa,-6,-1)=="000100",T,F) ## "000100" is water with no clouds
lsat<-lsat[lsat$qa_good==T,]
lsat<-lsat[is.na(lsat$B1)==F,]

##get location info from google earths .geo character string
x<-strsplit(lsat$".geo", "\\[")
x<-sapply(x, "[[", 2)
x<-strsplit(x, "\\]")
x<-sapply(x, "[[", 1)
x<-strsplit(x, ",")
lsat$long<-sapply(x, "[[", 1)
lsat$lat<-sapply(x, "[[", 2)
head(lsat)

##get landsat dates to feed to webtool
dates<-as.Date(unique(lsat$date))
##remove dates we already analyzed, if appending to a previous file
if (append.out==T){
  old$date<-as.Date(old$date)
  dates<-dates[!(dates %in% old$date)]
}

##get a lat long to give to webtool, when using profile = "2" it returns an 
##average over 1 degree lat and long, so we only need one location if the data are all from the same area
lat<-as.numeric(lsat$lat[1])
long<-as.numeric(lsat$long[1])

##a dataframe to hold the results, with the right number of rows pre allocated
out_atmos<-data.frame("date"=dates, "Iu"=rep(NA, length(dates)), "Id"=rep(NA, length(dates)), "trans"=rep(NA, length(dates)))

for (i in seq_along(dates)){
    x<- atmos(date=dates[i], lat=lat, long=long)
    out_atmos[as.Date(out_atmos$date)==as.Date(x$date),"trans"]<-x$trans    
    out_atmos[as.Date(out_atmos$date)==as.Date(x$date),"Id"]<-x$Id
    out_atmos[as.Date(out_atmos$date)==as.Date(x$date),"Iu"]<-x$Iu
}

##optionally add in the data from previous files before writing out, append.out option was set above
if(append.out==T){
  out_atmos<-rbind(old, out_atmos) }

##order the rows by date and write out
out_atmos<-arrange(out_atmos, date)
out_atmos<-out_atmos[complete.cases(out_atmos),]
write.csv(out_atmos, paste0(basedir,"output/atmospheric_correction_varsL8_gce.csv"), row.names=F)