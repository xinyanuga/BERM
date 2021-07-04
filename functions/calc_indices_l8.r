##calculate lsat 8 indices
##see https://landsat.usgs.gov/sites/default/files/documents/si_product_guide.pdf
##https://www.harrisgeospatial.com/docs/BroadbandGreenness.html
calc_index_l8<-function(lsat){
  ndvi<-(lsat$b5-lsat$b4)/(lsat$b5+lsat$b4) 
  indices<-data.frame(ndvi=ndvi)
  indices$ndmi<-(lsat$b1-lsat$b6)/(lsat$b1+lsat$b6) 
  indices$pheno<-(lsat$b5-lsat$b6)/(lsat$b5+lsat$b6) ##from TMII paper, is NDMI in the usgs link above
  indices$wdrvi <- (0.1*lsat$b5 - lsat$b4)/(0.1*lsat$b5 + lsat$b4)
  #EVI below is 2 b EVI and 3 b EVI
  indices$evi2 <- 2.5*(lsat$b5 - lsat$b4)/(lsat$b5 + 2.4*lsat$b4+1)
  indices$evi3 <- 2.5*(lsat$b5 - lsat$b4)/(lsat$b5 + 6*lsat$b4+7.5*lsat$b2+1)
  indices$savi <- 1.5*(lsat$b5 - lsat$b4)/(lsat$b5 + lsat$b4+0.5)
  indices$msavi <- (2*lsat$b5 + 1 - sqrt ((2*lsat$b5 + 1)^2 - 8*(lsat$b5 - lsat$b4)))/2
  indices$vari<- (lsat$b3 - lsat$b4)/(lsat$b3+lsat$b4+lsat$b2)
  ##green atomsphereically resistant index
  ##Gitelson, A., Y. Kaufman, and M. Merzylak. "Use of a Green Channel in Remote Sensing of Global Vegetation from EOS-MODIS." Remote Sensing of Environment 58 (1996): 289-298.
  indices$gari<- (lsat$b5 -(lsat$b3-1.7*(lsat$b2-lsat$b4)))/(lsat$b5 +(lsat$b3-1.7*(lsat$b2-lsat$b4)))
  ##green chlorophyll index
  ##Gitelson, A., Y. Gritz, and M. Merzlyak. "Relationships Between Leaf Chlorophyll Content and Spectral Reflectance and Algorithms for Non-Destructive Chlorophyll Assessment in Higher Plant Leaves." Journal of Plant Physiology 160 (2003): 271-282.
  indices$gci<-(lsat$b5/lsat$b3)-1
  ##alternative chlorophyll is derivation of ndvi, gndvi
  indices$gndvi<-(lsat$b5-lsat$b3)/(lsat$b5+lsat$b3) 
  ##index to calc LAI from EVI
  indices$lai<- 3.618*indices$evi3 -0.118
  ##nonlinear index-- linearizes relationship in ndvi
  ##Rouse, J., R. Haas, J. Schell, and D. Deering. Monitoring Vegetation Systems in the Great Plains with ERTS. Third ERTS Symposium, NASA (1973): 309-317.
  indices$nli<-(lsat$b5^2-lsat$b4)/(lsat$b5^2+lsat$b4) 
  indices$dvi<-lsat$b5-lsat$b4
  indices$rdvi<-(lsat$b5-lsat$b4)/sqrt(lsat$b5+lsat$b4)
  indices$sr<-(lsat$b5/lsat$b4)
  indices$mtvi1<-1.2*(1.2*(lsat$b5-lsat$b3)-2.5*(lsat$b4-lsat$b3))
  indices$mtvi2<-(1.5*(1.2*(lsat$b5-lsat$b3)-2.5*(lsat$b4-lsat$b3)))/sqrt((2*lsat$b5+1)^2-(6*lsat$b5-5*sqrt(lsat$b4))-0.5)
  
  return(indices)
}