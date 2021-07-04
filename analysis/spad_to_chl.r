
##script to estimate CHL concentration from SPAD readings

##read in the chlorophyll extraction data needed to calibrate the SPAD sensor readings and create easier column names
chl<-read.csv("data/seagrant_chl.csv", skip=10, header=F,  stringsAsFactors = F)
colnames(chl)<-c("labnum","site","plot","chla","chlb","chl")

##read in the raw vegetation plot level data, process the date field and create easier column names
veg<-read.csv("data/seagrant_vegplot.csv")
veg$date<-as.Date(veg$date, format= "%m-%d-%Y")
veg$spad<-veg$chlorophyll
##subset the veg data to the dates corresponding to the chlorophyll extraction analysis
veg<-veg[veg$date>"2016-10-01"&veg$date<"2016-11-01",]
##subset only to the veg columns we need for this analysis
veg<-veg[, c("plot", "spad")]

##combine the vegetation and chlorophyll extraction data
chl<-merge(chl, veg, by="plot")
## remove original runs for rerun data
chl<-chl[-c(17,39,41),]

summary(spad<-lm(chl~spad, data=chl))
#plot(chl$spad, chl$chl, ylab= "Chlorophyll (mg/g)", xlab="SPAD", cex.lab=1.5, pch=19)
#abline(spad)

#write.csv(chl, "output/chl_calibration.csv", row.names = F)

