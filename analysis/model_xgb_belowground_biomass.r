##after running the process groundtruth dataset for xgboost belowround biomass.r script, 
##run this script to actually build the model
##this model uses a nested resampling workflow with spatial crossvalidation to train, tune, build, and then test the model
##note that predicting novel pixels code (eg site level data) is at the end of this script

library(xgboost); library("paradox");  library("mlr3");library("mlr3learners"); 
library("data.table");  library("caret"); library("tidyverse"); library("mlr3tuning"); library("randomizr")

##output and save new models and data sets?? Set to True/False below. 
##Note that True will overwrite the existing models and output datasets, so you'll want a backup if you don't intend this
save_model<-FALSE
#save_model<-TRUE


data$id<-1:nrow(data)
##subset the data to dates where we have at least a year of antecedant conditions, 
##and remove dates without groundtruth data for this model
data<-data[data$date>"2014-06-01"&data$date<"2019-05-15",]

data$group<-paste(data$site, data$date)

##we'll be creating training and testing datasets as part of that process, 
##so let's ensure we save a copy of the data in memory by creating a new data set called train
train<-data
train$mo<-as.numeric(format(train$date, "%m"))  

##select to just the vars we want to consider as predictors for training the model
train<-dplyr::select(train, id,date,site, year, group,bgm2.core.stemscaled,
                     elevation, greendoy, growingday,
                     predag.allom.l8, predchl.l8, predlai.l8, photo, predn.l8,
                     lagg_ag_1, lagg_chl_1, lagg_lai_1, lagg_photo_1, lagg_N_1, lagg_ndvi_1,
                     roll_ag_2, roll_chl_2, roll_lai_2, roll_photo_2, roll_N_2, roll_ndvi_2,
                     roll_ag_3, roll_chl_3, roll_lai_3, roll_photo_3, roll_N_3, roll_ndvi_3,
                     roll_ag_4, roll_chl_4, roll_lai_4, roll_photo_4, roll_N_4, roll_ndvi_4,
                     roll_ag_5, roll_chl_5, roll_lai_5, roll_photo_5, roll_N_5, roll_ndvi_5,
                     diff_ag_1, diff_chl_1, diff_lai_1, diff_photo_1, diff_N_1, diff_ndvi_1,
                     diff_ag_2, diff_chl_2, diff_lai_2, diff_photo_2, diff_N_2, diff_ndvi_2,
                     diff_ag_3, diff_chl_3, diff_lai_3, diff_photo_3, diff_N_3, diff_ndvi_3,
                     diff_ag_4, diff_chl_4, diff_lai_4, diff_photo_4, diff_N_4, diff_ndvi_4,
                     diff_ag_5, diff_chl_5, diff_lai_5, diff_photo_5, diff_N_5, diff_ndvi_5,
                     diff_N_growing, diff_chl_growing, diff_ag_growing, diff_lai_growing,
                     growag, growchl, growlai, growphoto, growN,
                     local_hiwater, local_lowater, flood_time, 
                     lagg_lochi_1, roll_lochi_2, roll_lochi_3, roll_lochi_4,  roll_lochi_5,
                     lagg_loclo_1, roll_loclo_2, roll_loclo_3, roll_loclo_4,  roll_loclo_5,
                     lst, lagg_lst_1, roll_lst_2, roll_lst_3, roll_lst_4, roll_lst_5,
                     diff_lst_2,diff_lst_3,diff_lst_4,diff_lst_5,
                     dayl_mean, par_tot, prcp_mean, tmax_mean, tmin_mean,  vp_mean,
                     lagg_par_1, roll_par_2, roll_par_3, roll_par_4, roll_par_5,
                     lagg_prcp_1, roll_prcp_2, roll_prcp_3, roll_prcp_4, roll_prcp_5,
                     lagg_tmax_1, roll_tmax_2, roll_tmax_3,roll_tmax_4,roll_tmax_5,
                     lagg_tmin_1, roll_tmin_2, roll_tmin_3,roll_tmin_4, roll_tmin_5,
                     lagg_vp_1, roll_vp_2, roll_vp_3,roll_vp_4,roll_vp_5
)

##remove rows with missing data
train<-train[complete.cases(train),]
train$rowid<-as.numeric(rownames(train))

##feature names are the variables that will be predictor candidates, field data, dates, and row descriptors should be excluded
feature.names<-names(train)
feature.names<-feature.names[!(feature.names %in% c("date", "site","year","id", "rowid", "set","val","bgm2.core.stemscaled", 
                                                    "bgm2.allometric.rootshoot", "bgm2.allometric.rootgreenshoot",
                                                    "bgm2.allo.new.rootshoot", "bgm2.allo.new.rootgreenshoot", "mo",
                                                    "group"))]
feature.names

##Set up the resampling scheme for nested resampling
##outer is the number of outer crossvalidations
outer<-5
##note that for this workflow, inner models are just 1 model, 
##eg we divide the outer model training data just once into a single inner training and testing model 
##the inner training model is for hyperparameter tunning and feature selection, 
##which can be tested against an inner model testing set; 
##This preserves the outer testing set for final model fitting

##set up clusters and blocks, used to keep certain observations together during nested resampling
clusters <-train$group ##note, group is paste(site, date)
blocks <-  train$site
rep<-1:outer
vars<-paste0("cv",rep)
cvs<-matrix(data=NA,nrow=nrow(train),ncol=length(vars))
###create the cross valdidated folds, which will keep obsers from same site and date together####
set.seed(29382) 
for( i in 1: outer){
  cvs[,i]<- block_and_cluster_ra(blocks = blocks,
                                 clusters = clusters,
                                 prob_each = c(.35, .65))
}                          
cvs[1:4,]
cvs<-data.frame(cvs); names(cvs)<-vars
head(table(clusters, cvs$cv1));head(table(blocks, cvs$cv1))
head(table(clusters, cvs$cv2));head(table(blocks, cvs$cv2))
#take a look at the clusters; blocks
table(cvs$cv1)
table(blocks, cvs$cv1)
table(blocks, cvs$cv2)
table(blocks, cvs$cv3)

##create columns in the larger dataset that will be named "train.1" or "train.2" etc, 
##and that will have 0 or 1 if the row is part of that training set, where 0 means it's in the testing set; 
##this will provide a record for later access
for (i in 1:outer){
  train[[paste0("train.", i)]] <- cvs[,i]
}

for (i in 1:outer) {
  data[[paste0("train.", i)]]<-NA
  data[[paste0("train.", i)]] [data$id %in% train$id[train[[paste0("train.", i)]] ==1 ] ]<-1
  data[[paste0("train.", i)]] [data$id %in% train$id[train[[paste0("train.", i)]] ==0 ] ]<-0
}
table(data$train.1);table(data$train.2)
table(cvs$cv1);table(cvs$cv2)
train[train$site=="ugami"&train$train.2==0,"date"]
train[train$site=="fluxb"&train$train.2==0,"date"]
train[train$site=="fluxa"&train$train.2==1,"date"]

####Select the features to use####
###XGBoost####
resp<-"bgm2.core.stemscaled" ## the y variable

##create a cut off for feature importance across the corss validations; if the feature is less important it will be discarded
cut.off<-0.005
rep<-1:outer
vars<-paste0("train.",rep)
feats.all<-matrix(data=NA, nrow=length(feature.names), ncol=length(vars))
feats.all<-data.frame(feats.all); names(feats.all)<-vars
feats.all$features<-feature.names

##create a list of parameters to focus on tuning, must specify type (Int=Integer, Dbl=continuous number, Fct=factor/character)
tune_ps = ParamSet$new(list(
  #ParamFct$new("booster",  levels = c("gbtree","gblinear")),
  ParamInt$new("gamma", lower = 3, upper = 10),
  ParamDbl$new("eta", lower = 0.05, upper = 0.5),
  ParamInt$new("max_depth", lower = 4, upper = 11),
  ParamInt$new("min_child_weight", lower = 1, upper = 10),
  ParamDbl$new("subsample", lower = 0.25, upper = 0.8),
  ParamDbl$new("colsample_bytree", lower = 0.25, upper = 0.8)
))
tune_ps

####go through all CVs to get average features #####

for (i in 1:outer) {
  tr<-train [train[[paste0("train.", i)]] ==1, ]
  traintask = TaskRegr$new(id = "bgbiomass", backend = tr[,c(resp, "group", feature.names)], target = resp)
  traintask$select(feature.names)
  ##set up learner
  ##select a resampling stategy for the datasets and preformance measure to use to tune the paramters
  ##set up how the tuner should iterate through the hyperparamters to try; 
  ##set up how the tuning should stop: here it stops after 20 iterations
  ##see other options at: https://mlr3book.mlr-org.com/tuning.html
  learner = lrn("regr.xgboost", predict_type = "response",  booster="gbtree",  
                objective   = "reg:squarederror", nrounds=300, verbose=0)
  #resamp = rsmp("holdout")
  resamp=rsmp("holdout", ratio=0.85)
  measure = msr("regr.mse")
  tuner = mlr3tuning::tnr("grid_search", resolution = 20)
  #tuner = tnr("random_search")
  #terminator = trm("evals", n_evals = 10)
  terminator = trm("stagnation", iters=3, threshold=0.1)
  
  instance = mlr3tuning::TuningInstanceSingleCrit$new(
    task = traintask,
    learner = learner,
    resampling = resamp,
    measure = measure,
    search_space = tune_ps,
    terminator = terminator
  )
  instance
  tuner = mlr3tuning::tnr("grid_search", resolution = 20)
  #tuner = tnr("random_search")
  
  ##run the optimization
  tuner$optimize(instance)
  
  ##see the result, eg, the best configuration of parameters found during the iterations
  instance$result_learner_param_vals
  #see the benchmark results of all evaluations
  instance$archive$data()
  ##store the best parameters back in the learner
  learner$param_set$values = instance$result_learner_param_vals
  learner$train(traintask)
  
  ##gather the features to use
  y=learner$importance(); x=names(y)
  mar=c(3,15,2,2)
  features.all<-names(y[y>cut.off])
  features.all
  
  if(length(features.all)<length(feature.names)){
    feats.all[,i]<- c(features.all, rep(NA, length(feature.names)-length(features.all)))
  } else{
    feats.all[,i]<- features.all
  }
}

##summarize the features from the gathered features
features.all<-data.frame(var=feature.names, train.1=NA, train.2=NA, train.3=NA)
for( i in 1: outer){
  features.all [[paste0("train.", i)]]<-ifelse(features.all$var %in% feats.all[[paste0("train.", i)]], 1,0)
}

##create objects that record which features to retain and summarizes their importance across the crossvalidated trials
features.all$vote<-rowSums(features.all[,2:(outer+1)])
features.all.use<-features.all$var[features.all$vote>0.5*outer]

features<-features.all.use
features

#####Train the models with the selected features#####
###train model 1#####
resp<-"bgm2.core.stemscaled" ## the y variable
tr<-train[train$train.1==1,]
traintask = TaskRegr$new(id = "bgbiomass", backend = tr[,c(resp,  features)], target = resp)
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
instance$result_learner_param_vals
param <- instance$result_x_domain

resp<-tr$bgm2.core.stemscaled
clf <- xgboost(data        = data.matrix(tr[,features]),
               label       = resp,
               booster="gbtree",
               nrounds     = 125,
               params = param,
               verbose=0)

mat <- xgb.importance (feature_names = features,model = clf)
xgb.plot.importance (importance_matrix = mat[1:30]) 

##plot the training data predicted vs measured
p <-(as.integer(round(predict(clf, data.matrix(train[train$train.1==1,features])))))
plot(train$bgm2.core.stemscaled[train$train.1==1], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.1==1]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.1==1]))), col=factor(train$site[train$train.1==1])); abline(0,1)
sqrt((sum((train$bgm2.core.stemscaled[train$train.1==1]-p)^2, na.rm=T))/length(p))

##plot the testing data predicted vs measured
p <-(as.integer(round(predict(clf, data.matrix(train[train$train.1==0,features])))))
plot(train$bgm2.core.stemscaled[train$train.1==0], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.1==0]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.1==0]))), col=factor(train$site[train$train.1==0])); abline(0,1)
rmset<-sqrt((sum((train$bgm2.core.stemscaled[train$train.1==0]-p)^2, na.rm=T))/length(p))
print(rmset)

h<-train[train$train.1==0,]
h$p <-(as.integer(round(predict(clf, data.matrix(train[train$train.1==0,features])))))
cor(h$bgm2.core.stemscaled[h$site=="ugami"], h$p[h$site=="ugami"])
cor(h$bgm2.core.stemscaled[h$site=="fluxb"], h$p[h$site=="fluxb"])
cor(h$bgm2.core.stemscaled[h$site=="skida"], h$p[h$site=="skida"])

##save the model as it's own object
xgmod1<-clf
#clf<-xgmod1

#### Train set 2 #####
tr<-train [train$train.2==1, ]
resp<-"bgm2.core.stemscaled" ## the y variable
traintask = TaskRegr$new(id = "bgbiomass", backend = tr[,c(resp, features)], target = resp)
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
instance$result_learner_param_vals
param <- instance$result_x_domain


##Code with early stopping rounds based on hold out set performance
resp<-tr$bgm2.core.stemscaled
clf <- xgboost(data        = data.matrix(tr[,features]),
               label       = resp,
               booster="gbtree",
               nrounds     = 125,
               params = param, 
               verbose=0)

mat <- xgb.importance (feature_names = features,model = clf)
xgb.plot.importance (importance_matrix = mat[1:30]) 
p <-(as.integer(round(predict(clf, data.matrix(train[train$train.2==1,features])))))
plot(train$bgm2.core.stemscaled[train$train.2==1], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.2==1]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.2==1]))), col=factor(train$site[train$train.2==1])); abline(0,1)
sqrt((sum((train$bgm2.core.stemscaled[train$train.2==1]-p)^2, na.rm=T))/length(p))

p <-(as.integer(round(predict(clf, data.matrix(train[train$train.2==0,features])))))
plot(train$bgm2.core.stemscaled[train$train.2==0], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.2==0]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.2==0]))), col=factor(train$site[train$train.2==0])); abline(0,1)
rmset<-sqrt((sum((train$bgm2.core.stemscaled[train$train.2==0]-p)^2, na.rm=T))/length(p))
print(rmset)

h<-train[train$train.2==0,]
h$p <-(as.integer(round(predict(clf, data.matrix(train[train$train.2==0,features])))))
cor(h$bgm2.core.stemscaled[h$site=="ugami"], h$p[h$site=="ugami"])
cor(h$bgm2.core.stemscaled[h$site=="fluxb"], h$p[h$site=="fluxb"])
cor(h$bgm2.core.stemscaled[h$site=="skida"], h$p[h$site=="skida"])

xgmod2<-clf
#clf<-xgmod2

####Train set 3 ####
train$train.3[train$train.3==0&train$site=="ugami"][which.max(train$bgm2.core.stemscaled[train$train.3==0&train$site=="ugami"])]<-1
tr<-train [train$train.3==1, ]
resp<-"bgm2.core.stemscaled" ## the y variable
traintask = TaskRegr$new(id = "bgbiomass", backend = tr[,c(resp, features)], target = resp)
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
instance$result_learner_param_vals
param <- instance$result_x_domain


##Code with early stopping rounds based on hold out set performance
resp<-tr$bgm2.core.stemscaled
clf <- xgboost(data        = data.matrix(tr[,features]),
               label       = resp,
               booster="gbtree",
               nrounds     = 125,
               params = param, 
               verbose=0)

mat <- xgb.importance (feature_names = features,model = clf)
xgb.plot.importance (importance_matrix = mat[1:40]) 
p <-(as.integer(round(predict(clf, data.matrix(train[train$train.3==1,features])))))
plot(train$bgm2.core.stemscaled[train$train.3==1], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.3==1]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.3==1]))), col=factor(train$site[train$train.3==1])); abline(0,1)
sqrt((sum((train$bgm2.core.stemscaled[train$train.3==1]-p)^2, na.rm=T))/length(p))

p <-(as.integer(round(predict(clf, data.matrix(train[train$train.3==0,features])))))
plot(train$bgm2.core.stemscaled[train$train.3==0], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.3==0]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.3==0]))), col=factor(train$site[train$train.3==0])); abline(0,1)
rmset<-sqrt((sum((train$bgm2.core.stemscaled[train$train.3==0]-p)^2, na.rm=T))/length(p))
print(rmset)

h<-train[train$train.3==0,]
h$p <-(as.integer(round(predict(clf, data.matrix(train[train$train.3==0,features])))))
cor(h$bgm2.core.stemscaled[h$site=="ugami"], h$p[h$site=="ugami"])
cor(h$bgm2.core.stemscaled[h$site=="fluxb"], h$p[h$site=="fluxb"])
cor(h$bgm2.core.stemscaled[h$site=="skida"], h$p[h$site=="skida"])

xgmod3<-clf
#clf<-xgmod3

####Train set 4 ####
tr<-train [train$train.4==1, ]
resp<-"bgm2.core.stemscaled" ## the y variable
traintask = TaskRegr$new(id = "bgbiomass", backend = tr[,c(resp, features)], target = resp)
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
instance$result_learner_param_vals
param <- instance$result_x_domain


##Code with early stopping rounds based on hold out set performance
resp<-tr$bgm2.core.stemscaled
clf <- xgboost(data        = data.matrix(tr[,features]),
               label       = resp,
               booster="gbtree",
               nrounds     = 125,
               params = param, 
               verbose=0)

mat <- xgb.importance (feature_names = features,model = clf)
xgb.plot.importance (importance_matrix = mat[1:30]) 
p <-(as.integer(round(predict(clf, data.matrix(train[train$train.4==1,features])))))
plot(train$bgm2.core.stemscaled[train$train.4==1], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.4==1]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.4==1]))), col=factor(train$site[train$train.4==1])); abline(0,1)
sqrt((sum((train$bgm2.core.stemscaled[train$train.4==1]-p)^2, na.rm=T))/length(p))

p <-(as.integer(round(predict(clf, data.matrix(train[train$train.4==0,features])))))
plot(train$bgm2.core.stemscaled[train$train.4==0], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.4==0]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.4==0]))), col=factor(train$site[train$train.4==0])); abline(0,1)
rmset<-sqrt((sum((train$bgm2.core.stemscaled[train$train.4==0]-p)^2, na.rm=T))/length(p))
print(rmset)

h<-train[train$train.4==0,]
h$p <-(as.integer(round(predict(clf, data.matrix(train[train$train.4==0,features])))))
cor(h$bgm2.core.stemscaled[h$site=="ugami"], h$p[h$site=="ugami"])
cor(h$bgm2.core.stemscaled[h$site=="fluxb"], h$p[h$site=="fluxb"])
cor(h$bgm2.core.stemscaled[h$site=="skida"], h$p[h$site=="skida"])

xgmod4<-clf
#clf<-xgmod4

####Train set 5####
tr<-train [train$train.5==1, ]
resp<-"bgm2.core.stemscaled" ## the y variable
traintask = TaskRegr$new(id = "bgbiomass", backend = tr[,c(resp, features)], target = resp)
tuner$optimize(instance)

##see the result, eg, the best configuration of parameters found during the iterations
instance$result_learner_param_vals
param <- instance$result_x_domain


##Code with early stopping rounds based on hold out set performance
resp<-tr$bgm2.core.stemscaled
clf <- xgboost(data        = data.matrix(tr[,features]),
               label       = resp,
               booster="gbtree",
               nrounds     = 125,
               params = param, 
               verbose=0)

mat <- xgb.importance (feature_names = features,model = clf)
xgb.plot.importance (importance_matrix = mat[1:30]) 
p <-(as.integer(round(predict(clf, data.matrix(train[train$train.5==1,features])))))
plot(train$bgm2.core.stemscaled[train$train.5==1], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.5==1]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.5==1]))), col=factor(train$site[train$train.5==1])); abline(0,1)
sqrt((sum((train$bgm2.core.stemscaled[train$train.5==1]-p)^2, na.rm=T))/length(p))

p <-(as.integer(round(predict(clf, data.matrix(train[train$train.5==0,features])))))
plot(train$bgm2.core.stemscaled[train$train.5==0], p, ylim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.5==0]))),
     xlim=c(0,max(c(p, train$bgm2.core.stemscaled[train$train.5==0]))), col=factor(train$site[train$train.5==0])); abline(0,1)
rmset<-sqrt((sum((train$bgm2.core.stemscaled[train$train.5==0]-p)^2, na.rm=T))/length(p))
print(rmset)

h<-train[train$train.5==0,]
h$p <-(as.integer(round(predict(clf, data.matrix(train[train$train.5==0,features])))))
cor(h$bgm2.core.stemscaled[h$site=="ugami"], h$p[h$site=="ugami"])
cor(h$bgm2.core.stemscaled[h$site=="fluxb"], h$p[h$site=="fluxb"])
cor(h$bgm2.core.stemscaled[h$site=="skida"], h$p[h$site=="skida"])


xgmod5<-clf
#clf<-xgmod5

###gather features importance from the final fitted models from the outer model fitting####
rep<-1:outer
vars<-paste0("train.", rep)
feats.all.import<-matrix(data=NA, nrow=length(features), ncol=length(vars))
feats.all.import<-data.frame(feats.all.import); names(feats.all.import)<-vars
feats.all.import$features<-features

mat <- xgb.importance (feature_names = features,model = xgmod1)
x=mat$Feature; y=mat$Gain
import<-data.frame(features=x, import=as.numeric(y))
feats.all.import[match(import$features, feats.all.import$features),1]<- import$import

mat <- xgb.importance (feature_names = features,model = xgmod2)
x=mat$Feature; y=mat$Gain
import<-data.frame(features=x, import=as.numeric(y))
feats.all.import[match(import$features, feats.all.import$features),2]<- import$import

mat <- xgb.importance (feature_names = features,model = xgmod3)
x=mat$Feature; y=mat$Gain
import<-data.frame(features=x, import=as.numeric(y))
feats.all.import[match(import$features, feats.all.import$features),3]<- import$import

mat <- xgb.importance (feature_names = features,model = xgmod4)
x=mat$Feature; y=mat$Gain
import<-data.frame(features=x, import=as.numeric(y))
feats.all.import[match(import$features, feats.all.import$features),4]<- import$import

mat <- xgb.importance (feature_names = features,model = xgmod5)
x=mat$Feature; y=mat$Gain
import<-data.frame(features=x, import=as.numeric(y))
feats.all.import[match(import$features, feats.all.import$features),5]<- import$import

data$pred.1[!is.na(data$train.1)]<-(as.integer(round(predict(xgmod1, data.matrix(data[!is.na(data$train.1),features])))))
data$pred.2[!is.na(data$train.2)]<-(as.integer(round(predict(xgmod2, data.matrix(data[!is.na(data$train.2),features])))))
data$pred.3[!is.na(data$train.3)]<-(as.integer(round(predict(xgmod3, data.matrix(data[!is.na(data$train.3),features])))))
data$pred.4[!is.na(data$train.4)]<-(as.integer(round(predict(xgmod4, data.matrix(data[!is.na(data$train.4),features])))))
data$pred.5[!is.na(data$train.5)]<-(as.integer(round(predict(xgmod5, data.matrix(data[!is.na(data$train.5),features])))))

y<-dplyr::select(data, train.1, train.2, train.3,  train.4, train.5,
                 pred.1, pred.2, pred.3, pred.4, pred.5, bgm2.core.stemscaled)
rmse<-data.frame(set=1:(outer), train=rep(NA, (outer)), test=rep(NA, (outer)), 
                 trainmean=rep(NA, (outer)), trainmed=rep(NA, (outer)))
for (i in 1:outer){
  tr<-y[,c(i,i+outer,ncol(y))]
  colnames(tr)<-c("train", "pred", "resp")
  trainset<-tr[which(tr$train==1), ]
  testset<-tr[which(tr$train==0),]
  
  rmse$train[i]<-sqrt((sum((trainset$resp-trainset$pred)^2, na.rm=T))/length(trainset$resp[!is.na(trainset$resp)]))
  rmse$test[i]<-sqrt((sum((testset$resp-testset$pred)^2, na.rm=T))/length(testset$resp[!is.na(testset$resp)]))
  rmse$trainmean[i]<-mean(trainset$pred, na.rm=T)
  rmse$trainmed[i]<-median(trainset$pred, na.rm=T)
}
summary(rmse)
rmse<-arrange(rmse,set)
rmse

data$predbg.xgb<-rowMeans(data[,c("pred.1","pred.2", "pred.3", "pred.4", "pred.5")])
maxx<-max(c(data$bgm2.core.stemscaled, data$predbg.xgb), na.rm=T)
minn<-min(c(data$bgm2.core.stemscaled, data$predbg.xgb), na.rm=T)
sqrt((sum((data$bgm2.core.stemscaled-data$predbg.xgb)^2, na.rm=T))/length(data$predbg.xgb[!is.na(data$predbg.xgb)]))

jpeg(file="results/xgb_l8_bg_all_months.jpg", units="in",
     height=5.25, width=7.25, res=200)
par(mar=c(4,5,1,0.75), oma=c(0,0.2,0.5,0.5), mfrow=c(1,2))
layout(matrix(c(1,2),1,2,byrow = TRUE), c(1,2), TRUE)
boxplot(rmse[,   c("train","test")], 
        ylab =expression(paste("RMSE g ", m^-2, sep="")), 
        cex.lab=1.5, cex.axis=1.2)
mtext(side=3, adj=0.95, line=0.25, "A", cex=1.2)

plot( data$predbg.xgb, data$bgm2.core.stemscaled,
      xlab=expression(paste("predicted BG biomass g ", m^-2, sep="")), 
      ylab =expression(paste("measured BG biomass g ", m^-2, sep="")) , 
      ylim=c(minn,maxx),xlim=c(minn,maxx),  pch=c(19), 
      cex=c(1,0.8,1), col=trans.color[c(8,1,2,4)][factor(data$site)],
      cex.lab=1.5, cex.axis=1.2)
abline(0,1)
mtext(side=3, adj=0.05, line=-2, bquote(RMSE[train]: ~ .(round(mean(rmse$train, na.rm=T),1)) ~ g~m^-2), cex=1.2)
mtext(side=3, adj=0.05, line=-4, bquote(RMSE[test]: ~ .(round(mean(rmse$test, na.rm=T),1)) ~ g~m^-2), cex=1.2)
legend("bottomright", inset=c(0.0,0),legend=unique(data$site[!is.na(data$site)]), title=expression(underline(site)),title.adj = 0.25,
       col=trans.color[c(8,1,2,4)][factor(unique(data$site[!is.na(data$site)]))], pch=15, bty="n")
mtext(side=3, adj=0.97, line=0.25, "B", cex=1.2)
dev.off()

feat<-data.frame(all=c(features, rep(NA,length(feature.names)-length(features))), 
                 potential=feature.names)

##save the results
if(save_model==T){
  write_csv (data, "output/xgb_predicted_bgbiomass_seagrant_plots_landsat8.csv")
  write_csv (feat, "output/xgb_features.csv")
  write_csv (feats.all.import, "output/xgb_all_import_features.csv")
  write_csv (rmse, "output/xgb_model_rmse_seagrant_plots.csv")
  
  saveRDS(xgmod1, file = "output/xgboost_bgbiomass_1.rda")
  saveRDS(xgmod2, file = "output/xgboost_bgbiomass_2.rda")
  saveRDS(xgmod3, file = "output/xgboost_bgbiomass_3.rda")
  saveRDS(xgmod4, file = "output/xgboost_bgbiomass_4.rda")
  saveRDS(xgmod5, file = "output/xgboost_bgbiomass_5.rda")
  xgb.save(xgmod1,  "output/xgboost_bgbiomass_1")
  xgb.save(xgmod2,  "output/xgboost_bgbiomass_2")
  xgb.save(xgmod3,  "output/xgboost_bgbiomass_3")
  xgb.save(xgmod4,  "output/xgboost_bgbiomass_4")
  xgb.save(xgmod5,  "output/xgboost_bgbiomass_5")
  
}

##read in results previously saved, if needed
features<-read_csv("output/xgb_features.csv")
features<-features$all[!is.na(features$all)]
xgmod1<-xgb.load("output/xgboost_bgbiomass_1")
xgmod2<-xgb.load("output/xgboost_bgbiomass_2")
xgmod3<-xgb.load("output/xgboost_bgbiomass_3")
xgmod4<-xgb.load("output/xgboost_bgbiomass_4")
xgmod5<-xgb.load("output/xgboost_bgbiomass_5")

##Now we can use the models to predict novel data, typically at the site level
##make sure the novel data have run through the process_new_pixels_for_xgb.r script
newpixels<-read.csv("output/xgb_processed_seagrant_landsat8_allpixels.csv", stringsAsFactors = F)
newpixels$date<-as.Date(newpixels$date)

newout<-matrix(data=NA, nrow=nrow(newpixels), ncol=outer)
newout[,1]<-(as.integer(round(predict(xgmod1, data.matrix(newpixels[,features])))))
newout[,2]<-(as.integer(round(predict(xgmod2, data.matrix(newpixels[,features])))))
newout[,3]<-(as.integer(round(predict(xgmod3, data.matrix(newpixels[,features])))))
newout[,4]<-(as.integer(round(predict(xgmod4, data.matrix(newpixels[,features])))))
newout[,5]<-(as.integer(round(predict(xgmod5, data.matrix(newpixels[,features])))))

##save the predicted belowground biomass as the mean of the predictions from all 5 models
newpixels$predbg.xgb<-rowMeans(newout)
##set a floor for the predictions, though probably not even needed with XGB
newpixels$predbg.xgb<-ifelse(newpixels$predbg.xgb<0,50,newpixels$predbg.xgb)
newpixels$predag.allom.l8 <-ifelse(newpixels$predag.allom.l8<0,50,newpixels$predag.allom.l8)

##newpixels is too huge to keep all the variables, lets parse it down just to interesting ones for plotting results
newpixels<-dplyr::select(newpixels, pix,  year,  utm_east, utm_north, greenup, date, obs,  predn.l8,predchl.l8,
                         predlai.l8, predag.allom.l8, lst, ndvi, mo, greendoy,doy, growingday, photo, growag,  growN, 
                         diff_N_growing, diff_lai_growing, diff_chl_growing, diff_ag_growing,
                         growchl, growlai, growphoto, utm_east,utm_north,  elevation, water,hiwater, maxwater,  
                         lowater, local_hiwater, local_lowater,flood_time, dayl_mean,  prcp_mean,  srad_mean, par_sum,par_tot,
                         tmax_mean,  tmin_mean,  vp_mean, dayl_sum, prcp_sum,srad_sum,  tmax_sum,tmin_sum,vp_sum, 
                         predbg.xgb , flood_time)

##save the results
if(save_model==T){
  write_csv(newpixels, "output/xgb_predicted_bgbiomass_seagrant_landsat8_allpixels.csv")
}