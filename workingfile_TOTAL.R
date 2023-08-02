#Code used to generate forecasts of TOTAL
#For full methods and details see: Brodie et al., (in review)
#Data is deposited on dryad: (INSERT LINK)
#Steps below correspond to files stored on dryad


#------load librarys-------
library(ncdf4)
library(raster)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(glue)
library(zoo)
library(qmap)
library(ggpubr)
library(foreach)
library(measures)
library(lubridate)
library(verification)

#-----Step-1: FUNCTIONS: TPR, TNR, Accuracy------
#Skill Functions

#tp = true positive
#tn = true negative
#fp = false positive (forecast but not observed, i.e. predict a MHW but no MHW observed)
#fn = false negative (observed but not forecast, i.e. observed a MHW but didn't predict it)

tpr <- function(o,f,qo,qf){ #observed, forecast, observed quantile, and forecast quantile
  # to <- round(quantile(o,qo, na.rm=T),2) #observed threshold
  # tf <- round(quantile(f,qf, na.rm=T),2) #forecast threshold
  tp <- length(which(o>=qo & f>=qf))
  tn <- length(which(o<qo & f<qf))
  fp <- length(which(o<qo & f>qf))
  fn <- length(which(o>qo & f<qf))
  tpr <- tp/(tp+fn)
  return(tpr)
}

tnr <- function(o,f,qo,qf){
  # to <- round(quantile(o,qo, na.rm=T),2) #observed threshold
  # tf <- round(quantile(f,qf, na.rm=T),2) #forecast threshold
  tp <- length(which(o>=qo & f>=qf))
  tn <- length(which(o<qo & f<qf))
  fp <- length(which(o<qo & f>qf))
  fn <- length(which(o>qo & f<qf))
  tnr <- tn/(tn+fp)
  return(tnr)
}

acc <- function(o,f,qo,qf){
  # to <- round(quantile(o,qo, na.rm=T),2) #observed threshold
  # tf <- round(quantile(f,qf, na.rm=T),2) #forecast threshold
  tp <- length(which(o>=qo & f>=qf))
  tn <- length(which(o<qo & f<qf))
  fp <- length(which(o<qo & f>qf))
  fn <- length(which(o>qo & f<qf))
  acc <- (tp + tn) / sum(tp,tn,fp,fn)
  return(acc)
}





#-----Step 1a: create OBSERVED TOTAL: best case scenario-----
#Function
roms_clim_sst_function_full <- function(x,y){
  nc_hist <- nc_open(x)
  nc_nrt <- nc_open(y)
  #Get vars from netcdf
  lat <- ncvar_get(nc_hist,'lat')[1,]
  lon <- ncvar_get(nc_hist,'lon')[,1]
  year_hist <- ncvar_get(nc_hist,'year')
  month_hist <- ncvar_get(nc_hist,'month')
  year_nrt <- ncvar_get(nc_nrt,'year')
  month_nrt <- ncvar_get(nc_nrt,'month')
  #Find dims for TOTAl box: but use resolution of global models
  lat_dimS <- which(lat==round(30.8,0))
  lat_dimN <- which(lat==round(34.5,0))
  lon_dimW <- which(lon==round(-120.3,0))
  lon_dimE <- which(lon==-116)
  #Get sst
  dat_hist <- ncvar_get(nc_hist, 'sst')
  dat_nrt <- ncvar_get(nc_nrt, 'sst')
  
  #Get monthly SST in Box A for historical period (1981-2010) 
  obs_sst_hist <- as.data.frame(matrix(NA, 360, 5))
  colnames(obs_sst_hist) <- c("year","month","day","date","sst")
  obs_sst_hist$year <- c(year_hist)
  obs_sst_hist$month <- c(month_hist)
  obs_sst_hist$day <- "01"
  obs_sst_hist$date <- as.Date(paste0(obs_sst_hist$year,'-',obs_sst_hist$month,'-',obs_sst_hist$day), "%Y-%m-%d")
  counter = 1
  for (y in 1981:2010){
    for (m in 1:12){
      idx_hist <- which(year_hist==y & month_hist==m)
      sst_hist <- mean(dat_hist[lon_dimW:lon_dimE,lat_dimS:lat_dimN, idx_hist], na.rm=T)
      obs_sst_hist[counter,5] <- sst_hist
      counter=counter+1
    }
  }
  #Get monthly SST in Box A for NRT period (2011-2020) 
  obs_sst_nrt <- as.data.frame(matrix(NA, 123, 5))
  colnames(obs_sst_nrt) <- c("year","month","day","date","sst")
  obs_sst_nrt$year <- c(year_nrt)
  obs_sst_nrt$month <- c(month_nrt)
  obs_sst_nrt$day <- "01"
  obs_sst_nrt$date <- as.Date(paste0(obs_sst_nrt$year,'-',obs_sst_nrt$month,'-',obs_sst_nrt$day), "%Y-%m-%d")
  counter = 1
  for (y in 2011:2021){
    for (m in 1:12){
      idx_nrt <- which(year_nrt==y & month_nrt==m)
      sst_nrt <- mean(dat_nrt[lon_dimW:lon_dimE,lat_dimS:lat_dimN, idx_nrt], na.rm=T)
      obs_sst_nrt[counter,5] <- sst_nrt
      counter=counter+1
    }
  }
  #Merge temporal time periods, get clim and anomalies 
  obs_sst <- rbind(obs_sst_hist, obs_sst_nrt)
  clim <- as.data.frame(matrix(NA,nrow=12,ncol=2))
  colnames(clim) <- c("month", "sst_clim")
  counter=1
  for (m in 1:12){
    clims <- mean(obs_sst$sst[obs_sst$year>=1981 & obs_sst$year<=2020 & obs_sst$month==m], na.rm=T)
    clim[counter,1] <- m
    clim[counter,2] <- clims
    counter=counter+1
  }
  
  #Join monthly and clim
  obs_sst_all <- left_join(obs_sst,clim,by="month")
  obs_sst_all$obs_anom <- obs_sst_all$sst - obs_sst_all$sst_clim
  obs_sst_all$roll_obs_anom <- zoo::rollmean(obs_sst_all$obs_anom, k=6, fill=NA, align='right')
  return(obs_sst_all)
}
#Both historical and NRT data
obs_anom_nrt <- roms_clim_sst_function_full(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/wcra31_sst_monthly_1981_2010.nc',
                                            y = '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/wcnrt_sst_monthly_201101_202103.nc')
# obs_anom_nrt <- obs_anom_nrt[obs_anom_nrt$date<="2020-07-01",]
mean_anom_6month=obs_anom_nrt %>% filter(date=="2014-02-01"|date=="2014-03-01"|date=="2014-04-01"|
                                           date=="2014-05-01"|date=="2014-06-01"|date=="2014-07-01"|
                                           date=="2014-12-01"|date=="2015-01-01"|date=="2015-02-01"|
                                           date=="2015-03-01"|date=="2015-04-01"|date=="2015-05-01"|
                                           date=="2015-12-01"|date=="2016-01-01"|date=="2016-02-01"|
                                           date=="2016-03-01"|date=="2016-04-01"|date=="2016-05-01") %>%
  summarise(mean=min(obs_anom))  ## 6 months preceding closures

plot(obs_anom_nrt$date[ obs_anom_nrt$year<=2020],obs_anom_nrt$roll_obs_anom[ obs_anom_nrt$year<=2020], type='l')
quantile(obs_anom_nrt$roll_obs_anom[obs_anom_nrt$year<=2020], c(0.74), na.rm=T)
# abline(h=round(mean_anom_6month,2), col="red")
abline(h=round(0.41,2), col="red")
abline(h=round(0.40,2), col="red")
abline(h=round(0.77,2), col="red")

plot(obs_anom_nrt$date[ obs_anom_nrt$year<=2010],obs_anom_nrt$roll_obs_anom[ obs_anom_nrt$year<=2010], type='l')
quantile(obs_anom_nrt$roll_obs_anom[obs_anom_nrt$year<=2010], c(0.82), na.rm=T)
abline(h=round(0.41,2), col="red")

saveRDS(obs_anom_nrt, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1a/observed_anomalies_bestcase.rds')
saveRDS(round(mean_anom_6month,2), '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1a/observed_threshold_bestcase.rds')

#-----Step 1b: Create OBSERVED ROMS climatology: apples-to-apples scenario----
#Function
roms_clim_sst_function_hist <- function(x){
  nc <- nc_open(x)
  #Get vars from netcdf
  lat <- ncvar_get(nc,'lat')[1,]
  lon <- ncvar_get(nc,'lon')[,1]
  year <- ncvar_get(nc,'year')
  month <- ncvar_get(nc,'month')
  #Find dims for TOTAl box, but change to match global model resolution
  # Box: -120.3 - -116; 30.8-34.5
  lat_dimS <- which(lat==round(30.8,0))
  lat_dimN <- which(lat==round(34.5,0))
  lon_dimW <- which(lon==round(-120.3,0))
  lon_dimE <- which(lon==-116)
  #Get sst
  dat <- ncvar_get(nc, 'sst')
  
  #Get monthly climatology (1991-2010), value, & anomaly from 1981-2010 over box
  obs_sst_metrics <- as.data.frame(matrix(NA, 360, 7))
  colnames(obs_sst_metrics) <- c("year","month","day","date","sst","sst_clim","obs_anom")
  obs_sst_metrics$year <- year
  obs_sst_metrics$month <- month
  obs_sst_metrics$day <- "01"
  obs_sst_metrics$date <- as.Date(paste0(obs_sst_metrics$year,'-',obs_sst_metrics$month,'-',obs_sst_metrics$day), "%Y-%m-%d")
  counter = 1
  for (y in 1981:2010){
    for (m in 1:12){
      
      clim_idx <- which(obs_sst_metrics$year>=1991 & obs_sst_metrics$month==m)
      sst_idx <- which(obs_sst_metrics$year==y & obs_sst_metrics$month==m)
      
      clim <- mean(dat[lon_dimW:lon_dimE,lat_dimS:lat_dimN,clim_idx], na.rm=T) 
      sst <-  mean(dat[lon_dimW:lon_dimE,lat_dimS:lat_dimN,sst_idx], na.rm=T)
      
      obs_sst_metrics[counter,5] <- sst
      obs_sst_metrics[counter,6] <- clim
      obs_sst_metrics[counter,7] <- sst - clim
      counter=counter+1
    }
  }
  return(obs_sst_metrics)
}

#historical data
obs_anom <- roms_clim_sst_function_hist(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/wcra31_sst_monthly_1981_2010.nc')
obs_anom$roll_obs_anom <- zoo::rollmean(obs_anom$obs_anom, k=6, fill=NA, align='right') #6month rolling mean
plot(obs_anom$date,obs_anom$roll_obs_anom, type='l')
qa <- quantile(obs_anom_nrt$roll_obs_anom, probs=c(0.834), na.rm=T)
abline(h=round(qa,2), col="red")

saveRDS(obs_anom, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1b/observed_anomalies_apples.rds')
saveRDS(round(qa,2), '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1b/observed_threshold_apples.rds')

#-----Step 2a: Calculate forecast SSTA on Global models: best case scenario-----
#Loop through each global models file and get anomaly SSTA within TOTAL box
#Do this for each global model, member, lead time, month
# x = file name directory
# d = date of first forecast
sst_box_function_GlobalModels <- function(x,d){
  sst_box <- as.data.frame(matrix(NA, nrow=11400, ncol=6))
  colnames(sst_box) <- c("glob_model","ensemble","init_month","lead_month","forecast_date","sst_mean")
  sst_box$forecast_date <- as.Date(sst_box$forecast_date,"%Y-%m-%d")
  sst_box$init_month <- as.Date(sst_box$init_month,"%Y-%m-%d")
  counter=1
  for (f in 1:length(x)){
    nc <- x
    nc <- nc_open(nc)
    # print(nc)
    #Get vars from netcdf
    lat <- ncvar_get(nc,'lat')
    lon <- ncvar_get(nc,'lon')
    lead <- ncvar_get(nc,'lead')
    time <- ncvar_get(nc,'time') #months since 1960-01-01. I don't know how to convert this, so next few lines are a workaround
    t1 <- ymd(d)
    time <- time-time[1]
    init_time <- t1 %m+% months(time)
    clim_idx <- which(init_time>="1991-01-01" & init_time<="2020-12-31")
    clim_start <- clim_idx[1]
    clim_end <- tail(clim_idx,n=1)
    member <- ncvar_get(nc,'member')
    t_start <- which(init_time>="1981-01-01")[1]
    #Find dims for TOTAl box
    lat_dimS <- which(lat==round(30.8,0))
    lat_dimN <- which(lat==round(34.5,0))
    lon_dimW <- which(lon==round(239.7,0)) #-120.3+360=239.7
    lon_dimE <- which(lon==244) #-116+360=244
    #Get monthly climatology & make SSTA from 2003-2010, for each init month and each lead time & member
    dat <- ncvar_get(nc, 'sst')
    
    for (i in t_start:clim_end){ #years/month index from start to Dec-2010
      for (l in 1:length(lead)){ #index
        for (m in 1:length(member)){ #ensemble member index
          
          dat_monthly <- mean(dat[lon_dimW:lon_dimE,lat_dimS:lat_dimN,m,l,i], na.rm=T)
          
          #Get forecast month
          lead.i <- lead[l] - 0.5
          forecast_date <- ymd(init_time[i]) %m+% months(lead.i)
          
          #Write out data
          sst_box[counter,1] <- unlist(strsplit(x[f],"_"))[7] #global model
          sst_box[counter,2] <- m #ensemble
          sst_box[counter,3] <- init_time[i] #initialisation month index
          sst_box[counter,4] <- lead[l] #lead time 
          sst_box[counter,5] <- forecast_date #get forecast date
          sst_box[counter,6] <- dat_monthly
          counter=counter+1
        }
      }
    }
  }
  return(sst_box)
}

#Run  for global NNME #takes ~34mins
t1 <- Sys.time()
sst_box_CanCM4 <- sst_box_function_GlobalModels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_CanCM4i.nc', d= "1981-01-01") #10 ens
sst_box_COLA <- sst_box_function_GlobalModels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_COLA-RSMAS-CCSM4.nc', d= "1982-01-01")#10 ens
sst_box_NEMO <- sst_box_function_GlobalModels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_GEM-NEMO.nc', d= "1981-01-01")#10 ens
sst_box_GFDL <- sst_box_function_GlobalModels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_GFDL-SPEAR.nc', d= "1991-01-01")#15 ens
sst_box_NASA <- sst_box_function_GlobalModels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_NASA-GEOSS2S.nc', d= "1981-02-01") #4 ens
sst_box_NCEP <- sst_box_function_GlobalModels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_NCEP-CFSv2.nc', d= "1982-01-01") #24 ens
t2 <- Sys.time()
t2-t1
global_fcasts <- rbind(sst_box_CanCM4,sst_box_COLA,sst_box_NEMO,
                       sst_box_GFDL,sst_box_NASA,sst_box_NCEP)

#Generate monthly climatology
global_fcasts <- na.omit(global_fcasts)
global_fcasts$month <- lubridate::month(global_fcasts$forecast_date)
global_fcasts_clim <- global_fcasts[global_fcasts$forecast_date>="1981-01-01",]
clims <- global_fcasts_clim %>% group_by(glob_model, ensemble, lead_month, month) %>% summarise_at(vars("sst_mean"),funs(mean(., na.rm=T))) 
colnames(clims)[5] <- "sst_clim"
global_fcasts <- left_join(global_fcasts,clims[c(1,2,3,4,5)],by=c("glob_model","ensemble","month", "lead_month"))
global_fcasts$sst_anom <- global_fcasts$sst_mean - global_fcasts$sst_clim 
global_fcasts <- global_fcasts[global_fcasts$sst_mean>10,]

saveRDS(global_fcasts,"~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2a/global_fcasts.rds")

#Plot time-series of SST
temp <- global_fcasts
# temp2 <- temp
ggplot(data=temp, aes(x=forecast_date, y=sst_mean, group=glob_model))+
  geom_line(aes(col=glob_model))+
  theme(legend.position = "bottom")+
  ggtitle("Global SST forecasts")+
  facet_wrap(~lead_month)
ggplot(data=temp, aes(x=forecast_date, y=sst_clim, group=glob_model))+
  geom_line(aes(col=glob_model))+
  theme(legend.position = "bottom")+
  ggtitle("Global SST climatology (1991-2020) for global forecasts")+
  facet_wrap(~lead_month)
ggplot(data=temp, aes(x=init_month, y=sst_anom, group=glob_model))+
  geom_line(aes(col=glob_model))+
  theme(legend.position = "bottom")+
  ggtitle("Global SST_A forecasts")+
  facet_wrap(~lead_month)

##get SSTA thresholds for each lead time 
#First need to ensemble
ens <- global_fcasts %>% group_by(lead_month, init_month, forecast_date) %>% summarise_at("sst_anom",mean, na.rm=T)
ens_roll <- ens %>% group_by(lead_month) %>% summarise_at(vars("forecast_date","init_month","sst_anom"),funs(rollmean(.,k=6,fill=NA,align='right')))
ens_thresholds <- ens_roll %>% group_by(lead_month) %>%  summarise(threshold=quantile(sst_anom, c(0.74),na.rm=T))
 
#Get 6-month rolling mean 
global_ensemble_rolling <- as.data.frame(matrix(NA, nrow=100,ncol=3))
colnames(global_ensemble_rolling) <- c("init_month","lead","sst_anom_rolling")
global_ensemble_rolling$init_month <- as.Date(global_ensemble_rolling$init_month)
counter=1
dates <- as.Date(unique(ens$init_month),"%Y-%m-$d")
for (d in 1:length(dates)){
  for (l in 0.5:5.5){
    dt <- dates[d]
    sst_anom_fcast_rolling <- mean(ens$sst_anom[ens$init_month==dt &
                                                  ens$lead_month>=l & ens$lead_month<=l+6],na.rm=T)
    global_ensemble_rolling[counter,1] <- as.Date(dt)
    global_ensemble_rolling[counter,2] <- l
    global_ensemble_rolling[counter,3] <- sst_anom_fcast_rolling
    counter=counter+1
  }
}
global_ensemble_rolling
global_ensemble_rolling$month <- lubridate::month(global_ensemble_rolling$init_month)
#Get Quantiles
quants <- global_ensemble_rolling %>% group_by(lead) %>% summarise_at("sst_anom_rolling", funs(quantile(.,c(0.74),na.rm=T)))
colnames(quants)[2] <- "ssta_quantile"
quants$ssta_quantile <- round(quants$ssta_quantile,2)

global_ensemble_rolling$forecast_date <- global_ensemble_rolling$init_month %m+%  months(global_ensemble_rolling$lead+0.5)
dt <- left_join(global_ensemble_rolling, quants, "lead")
ggplot(dt, aes(x=forecast_date, y=sst_anom_rolling))+
  geom_line()+
  geom_hline(aes(yintercept = ssta_quantile), col="red")+
  facet_wrap(~lead)

saveRDS(quants, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2a/thresholds_globalforecasts_bestcase.rds')



#-----Step 2b: Calculate forecast SSTA on Global models: apples-to-apples scenario -----

global_fcasts <- readRDS("~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2a/global_fcasts.rds")
global_fcasts <- global_fcasts[,c(1:7)]
global_fcasts_clim <- global_fcasts[global_fcasts$forecast_date>="1981-01-01" & global_fcasts$forecast_date<"2011-01-01",]
clims <- global_fcasts_clim %>% group_by(glob_model, ensemble, lead_month, month) %>% summarise_at(vars("sst_mean"),funs(mean(., na.rm=T))) 
colnames(clims)[5] <- "sst_clim"
global_fcasts <- left_join(global_fcasts,clims[c(1,2,3,4,5)],by=c("glob_model","ensemble","month", "lead_month"))
global_fcasts$sst_anom <- global_fcasts$sst_mean - global_fcasts$sst_clim 
global_fcasts <- global_fcasts[global_fcasts$sst_mean>10,]
global_fcasts_apples <- global_fcasts[global_fcasts$forecast_date<="2010-12-01",]
saveRDS(global_fcasts_apples,"~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2b/global_fcasts_apples.rds")

##get SSTA thresholds for each lead time 
#First need to ensemble
ens <- global_fcasts_apples %>% group_by(lead_month, init_month, forecast_date) %>% summarise_at("sst_anom",mean, na.rm=T)

#Get 6-month rolling mean 
global_ensemble_rolling <- as.data.frame(matrix(NA, nrow=100,ncol=3))
colnames(global_ensemble_rolling) <- c("init_month","lead","sst_anom_rolling")
global_ensemble_rolling$init_month <- as.Date(global_ensemble_rolling$init_month)
counter=1
dates <- as.Date(unique(ens$init_month),"%Y-%m-$d")
for (d in 1:length(dates)){
  for (l in 0.5:5.5){
    dt <- dates[d]
    sst_anom_fcast_rolling <- mean(ens$sst_anom[ens$init_month==dt &
                                                  ens$lead_month>=l & ens$lead_month<=l+6],na.rm=T)
    global_ensemble_rolling[counter,1] <- as.Date(dt)
    global_ensemble_rolling[counter,2] <- l
    global_ensemble_rolling[counter,3] <- sst_anom_fcast_rolling
    counter=counter+1
  }
}
global_ensemble_rolling
global_ensemble_rolling$month <- lubridate::month(global_ensemble_rolling$init_month)
#Get Quantiles
quants <- global_ensemble_rolling %>% group_by(lead) %>% summarise_at("sst_anom_rolling", funs(quantile(.,c(0.82),na.rm=T)))
colnames(quants)[2] <- "ssta_quantile"
quants$ssta_quantile <- round(quants$ssta_quantile,2)

global_ensemble_rolling$forecast_date <- global_ensemble_rolling$init_month %m+%  months(global_ensemble_rolling$lead+0.5)
ggplot(global_ensemble_rolling, aes(x=forecast_date, y=sst_anom_rolling))+
  geom_line()+
  geom_hline(yintercept = 0.11, col="red")+
  facet_wrap(~lead)

saveRDS(quants, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2b/thresholds_globalforecasts_apples.rds')

#-----Step 2c: Calculate forecast SSTA on Global models: only-3 scenario -----

global_fcasts <- readRDS("~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2b/global_fcasts_apples.rds")
global_fcasts_only3 <- global_fcasts_apples[global_fcasts_apples$glob_model=="CanCM4i.nc",]
global_fcasts_only3 <- global_fcasts_only3[global_fcasts_only3$ensemble==2 | global_fcasts_only3$ensemble==8 | global_fcasts_only3$ensemble==10,]

saveRDS(global_fcasts_only3,"~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2c/global_fcasts_only3.rds")

##get SSTA thresholds for each lead time 
#First need to ensemble
ens <- global_fcasts_only3 %>% group_by(lead_month, init_month, forecast_date) %>% summarise_at("sst_anom",mean, na.rm=T)

#Get 6-month rolling mean 
global_ensemble_rolling <- as.data.frame(matrix(NA, nrow=100,ncol=3))
colnames(global_ensemble_rolling) <- c("init_month","lead","sst_anom_rolling")
global_ensemble_rolling$init_month <- as.Date(global_ensemble_rolling$init_month)
counter=1
dates <- as.Date(unique(ens$init_month),"%Y-%m-$d")
for (d in 1:length(dates)){
  for (l in 0.5:5.5){
    dt <- dates[d]
    sst_anom_fcast_rolling <- mean(ens$sst_anom[ens$init_month==dt &
                                                  ens$lead_month>=l & ens$lead_month<=l+6],na.rm=T)
    global_ensemble_rolling[counter,1] <- as.Date(dt)
    global_ensemble_rolling[counter,2] <- l
    global_ensemble_rolling[counter,3] <- sst_anom_fcast_rolling
    counter=counter+1
  }
}
global_ensemble_rolling
global_ensemble_rolling$month <- lubridate::month(global_ensemble_rolling$init_month)
#Get Quantiles
quants <- global_ensemble_rolling %>% group_by(lead) %>% summarise_at("sst_anom_rolling", funs(quantile(.,c(0.82),na.rm=T)))
colnames(quants)[2] <- "ssta_quantile"
quants$ssta_quantile <- round(quants$ssta_quantile,2)

global_ensemble_rolling$forecast_date <- global_ensemble_rolling$init_month %m+%  months(global_ensemble_rolling$lead+0.5)
ggplot(global_ensemble_rolling, aes(x=forecast_date, y=sst_anom_rolling))+
  geom_line()+
  geom_hline(yintercept = 0.11, col="red")+
  facet_wrap(~lead)

saveRDS(quants, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2c/thresholds_globalforecasts_only3.rds')

#-----Step 3b: Calculate forecast SSTA on Downscaled models: apples-to-apple scenario-----
#Loop through each ens file and get clim, mean, and anomaly SST within TOTAL box:
#What the file needs to contain: all years (29), all ensembles + mean (n=4), all lead times (n=12)

sst_box_function <- function(x){
  sst_box <- as.data.frame(matrix(NA, nrow=1392, ncol=5))
  colnames(sst_box) <- c("ensemble","init_month","lead_month","forecast_date","sst_mean")
  sst_box$forecast_date <- as.Date(sst_box$forecast_date,"%Y-%m-%d")
  sst_box$init_month <- as.Date(sst_box$init_month,"%Y-%m-%d")
  counter=1
  for (f in 1:length(x)){
    nc <- x[f]
    nc <- nc_open(nc)
    # print(nc)
    #Get vars from netcdf
    lat <- ncvar_get(nc,'lat')[1,]
    lon <- ncvar_get(nc,'lon')[,1]
    lead <- ncvar_get(nc,'lead_time')
    init_time <- as.Date(ncvar_get(nc,'init_time'), "1970-01-01")
    #Find dims for TOTAl box;     #Change to match global model domain
    lat_dimS <- which(lat==round(30.8,0))
    lat_dimN <- which(lat==round(34.5,0))
    lon_dimW <- which(lon==round(-120.3,0)) 
    lon_dimE <- which(lon==-116) 
    #Get data
    dat <- ncvar_get(nc, 'sst')
    
    for (l in 1:12){ #indices
      for (y in 1:29){ #indices
        #Clim has to be calculated for each lead time. Don't do it here, do it outside of function
        dat_annual <- mean(dat[lon_dimW:lon_dimE,lat_dimS:lat_dimN,l,y], na.rm=T)
        #Get initialisation date (1 or 7)
        init <- unlist(strsplit(x[f],"_"))[5]
        init_month <- ifelse(init=="init1",1,7)
        #Get forcast month
        lead.i <- lead[l] - 0.5
        forecast_date <- ymd(init_time[y]) %m+% months(lead.i)
        #Write out data
        sst_box[counter,1] <- unlist(strsplit(x[f],"_"))[4] #ensemble
        sst_box[counter,2] <- init_time[y]
        sst_box[counter,3] <- lead[l]  #lead month
        sst_box[counter,4] <- forecast_date #get forecast date
        sst_box[counter,5] <- dat_annual
        counter=counter+1
      }
    }
  }
  return(sst_box)
}

#Load files
init1_files <- list.files('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/ForecastFields/', full.names = T, pattern="init1")
init7_files <- list.files('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/ForecastFields/', full.names = T, pattern="init7")
sst_files_init1 <- init1_files[grep('sst',init1_files)]
sst_files_init7 <- init7_files[grep('sst',init7_files)]

#Run function
sst_box_init1 <- sst_box_function(x = sst_files_init1)
sst_box_init7 <- sst_box_function(x = sst_files_init7)
downscaled_fcasts <- rbind(sst_box_init1,sst_box_init7)
downscaled_fcasts <- downscaled_fcasts[downscaled_fcasts$ensemble!="ensmean",] #get rid of mean and calculate it later from SSTA

#Generate monthly climatology for each lead time
downscaled_fcasts <- na.omit(downscaled_fcasts)
downscaled_fcasts$month <- lubridate::month(downscaled_fcasts$forecast_date)
downscaled_fcasts_trim <- downscaled_fcasts[downscaled_fcasts$init_month>="1981-01-01",]
clims <- downscaled_fcasts %>% group_by(ensemble, lead_month, month) %>% summarise_at(vars("sst_mean"),funs(mean(., na.rm=T))) 
colnames(clims)[4] <- "sst_clim"
downscaled_fcasts <- left_join(downscaled_fcasts,clims,by=c("ensemble","month", "lead_month"))
downscaled_fcasts$sst_anom <- downscaled_fcasts$sst_mean - downscaled_fcasts$sst_clim 

ggplot(data=downscaled_fcasts, aes(x=forecast_date, y=sst_mean, group=ensemble))+
  geom_line(aes(col=ensemble))+
  theme(legend.position = "bottom")+
  ggtitle("Downscaled SST forecasts")+
  facet_wrap(~lead_month)
ggplot(data=downscaled_fcasts, aes(x=forecast_date, y=sst_clim, group=ensemble))+
  geom_line(aes(col=ensemble))+
  theme(legend.position = "bottom")+
  ggtitle("Downscaled SST climatology (1991-2010) for global forecasts")+
  facet_wrap(~lead_month)
ggplot(data=downscaled_fcasts, aes(x=forecast_date, y=sst_anom, group=ensemble))+
  geom_line(aes(col=ensemble))+
  theme(legend.position = "bottom")+
  ggtitle("Downscaled SST_A forecasts")+
  facet_wrap(~lead_month)


saveRDS(downscaled_fcasts, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_3b/downscaled_forecasts.rds')

#Make ensemble mean & find thresholds
downscaled_ensemble <-  downscaled_fcasts %>% group_by(init_month,lead_month,forecast_date,month) %>%
  summarise_at(vars("sst_mean","sst_clim", "sst_anom"),funs(mean(., na.rm=T)))

#Get 6-month rolling mean 
downscaled_ensemble_rolling <- as.data.frame(matrix(NA, nrow=100,ncol=3))
colnames(downscaled_ensemble_rolling) <- c("init_month","lead","sst_anom_rolling")
downscaled_ensemble_rolling$init_month <- as.Date(downscaled_ensemble_rolling$init_month)
counter=1
dates <- as.Date(unique(downscaled_ensemble$init_month),"%Y-%m-$d")
for (d in 1:length(dates)){
  for (l in 0.5:5.5){
    dt <- dates[d]
    sst_anom_fcast_rolling <- mean(downscaled_ensemble$sst_anom[downscaled_ensemble$init_month==dt &
                                                                  downscaled_ensemble$lead_month>=l & downscaled_ensemble$lead_month<=l+6],na.rm=T)
    downscaled_ensemble_rolling[counter,1] <- as.Date(dt)
    downscaled_ensemble_rolling[counter,2] <- l
    downscaled_ensemble_rolling[counter,3] <- sst_anom_fcast_rolling
    counter=counter+1
  }
}
downscaled_ensemble_rolling
downscaled_ensemble_rolling$month <- lubridate::month(downscaled_ensemble_rolling$init_month)
#Get Quantiles
quants <- downscaled_ensemble_rolling %>% group_by(lead) %>% summarise_at("sst_anom_rolling", funs(quantile(.,c(0.82),na.rm=T)))
colnames(quants)[2] <- "ssta_quantile"
quants$ssta_quantile <- round(quants$ssta_quantile,2)
saveRDS(quants, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_3b/ssta_total_forecast_quantiles_downscaled.rds')
#merge and plot
downscaled_ensemble_rolling_withQuantiles <- left_join(downscaled_ensemble_rolling,quants,"lead")
ggplot(data=downscaled_ensemble_rolling_withQuantiles)+
  geom_line(aes(x=init_month, y=sst_anom_rolling))+
  geom_hline(aes(yintercept=ssta_quantile),col="red")+
  facet_wrap(~lead)


#-----Step 4a: Calculate SSTA for 6-month period: best case scenario-------
#load in data & join with obs
obs_anom <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1a/observed_anomalies_bestcase.rds')
glob <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2a/global_fcasts.rds')

#prepare obs for merge
obs_anom$sst_obs <- obs_anom$sst
obs_anom$sst_clim_obs <- obs_anom$sst_clim
obs_anom$sst_anom_obs <- obs_anom$obs_anom
obs_anom$sst_anom_obs_rolling <- obs_anom$roll_obs_anom
obs_anom$forecast_date <- obs_anom$date
obs_anom <- obs_anom[,c(9,10,11,12,13)]

all_data <- left_join(glob,obs_anom,"forecast_date")
all_data$year <- lubridate::year(all_data$forecast_date)
all_data <- na.omit(all_data)

globs <- c("CanCM4i.nc","COLA-RSMAS-CCSM4.nc", "GEM-NEMO.nc",
            "GFDL-SPEAR.nc", "NASA-GEOSS2S.nc", "NCEP-CFSv2.nc")
df_metric <- as.data.frame(matrix(NA, nrow=63, ncol=10))
colnames(df_metric) <- c("glob_model","ensemble","lead", "year", "june_fcast", "june_obs", "jul_fcast", "jul_obs", "aug_fcast", "aug_obs")
counter=1
t1 <- Sys.time() #Takes 12 mins
for (g in globs){
  for (l in 0.5:5.5){ #start lead time, which then adds 6 months to get lead end time. 
    for (m in 1:24){ #1:24
      for (y in 1981:2020){  
        print(glue("model: {g}; lead {l}; ens {m}; year {y}"))
        #Get mean over Dec-may
        # get forecast dates for Jun: dec- may
        dec <- as.Date(glue('{y}-12-01'))
        may <- as.Date(glue('{y+1}-05-01')) #do it this way to make sure you get correct year
        init_date <-  ymd(glue('{y}-12-01')) %m-% months(l-0.5)
        
        jun_fcast <- mean(all_data$sst_anom[all_data$glob_model==g &
                                              all_data$ensemble==m &
                                              all_data$init_month==init_date&
                                              all_data$forecast_date>=dec & 
                                              all_data$forecast_date<=may], na.rm=T)
        jun_obs <- mean(all_data$sst_anom_obs[ all_data$glob_model==g &
                                                 all_data$ensemble==m &
                                                 all_data$init_month==init_date&
                                                 all_data$forecast_date>=dec & 
                                                 all_data$forecast_date<=may], na.rm=T)
        
        #Get mean over Jan-Jun
        init_date <-  ymd(glue('{y+1}-01-01')) %m-% months(l-0.5)
        jul_fcast <- mean(all_data$sst_anom[all_data$glob_model==g &
                                              all_data$ensemble==m &
                                              all_data$init_month==init_date&
                                              all_data$month<=6], na.rm=T)
        jul_obs <- mean(all_data$sst_anom_obs[all_data$glob_model==g &
                                                        all_data$ensemble==m &
                                                        all_data$init_month==init_date&
                                                        all_data$month<=6], na.rm=T)
        
        #Get mean over Feb-Jul
        init_date <-  ymd(glue('{y+1}-02-01')) %m-% months(l-0.5)
        aug_fcast <- mean(all_data$sst_anom[all_data$glob_model==g &
                                              all_data$ensemble==m &
                                              all_data$init_month==init_date&
                                              all_data$month<=7 & all_data$month>=2],na.rm=T)
        aug_obs <- mean(all_data$sst_anom_obs[all_data$glob_model==g &
                                                all_data$ensemble==m &
                                                all_data$init_month==init_date&
                                                all_data$month<=7 & all_data$month>=2],na.rm=T)
        #write out data
        df_metric[counter,1] <- g
        df_metric[counter,2] <- m
        df_metric[counter,3] <- l
        df_metric[counter,4] <- y+1
        df_metric[counter,5] <- jun_fcast
        df_metric[counter,6] <- jun_obs
        df_metric[counter,7] <- jul_fcast
        df_metric[counter,8] <- jul_obs
        df_metric[counter,9] <- aug_fcast
        df_metric[counter,10] <- aug_obs
        counter=counter+1
      }
    }
  }
}
t2 <- Sys.time()
t2-t1
df_metric

saveRDS(df_metric,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_4a/total_metric_globalforecasts_1980-2020_bestcase.rds')

ggplot(data=df_metric,aes(x=june_fcast, y=june_obs))+
  geom_point(aes(col=glob_model),alpha=0.6)+
  ggtitle("June Rule")+
  stat_cor(method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0,col="dark grey")+
  geom_vline(xintercept=0, col="dark grey")+
  facet_wrap(~lead)

ggplot(data=df_metric,aes(x=jul_fcast, y=jul_obs))+
  geom_point(aes(col=glob_model),alpha=0.6)+
  ggtitle("July Rule")+
  stat_cor(method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0,col="dark grey")+
  geom_vline(xintercept=0, col="dark grey")+
  facet_wrap(~lead)

ggplot(data=df_metric,aes(x=aug_fcast, y=aug_obs))+
  geom_point(aes(col=glob_model),alpha=0.6)+
  ggtitle("August Rule")+
  stat_cor(method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0,col="dark grey")+
  geom_vline(xintercept=0, col="dark grey")+
  facet_wrap(~lead)

#Another way of looking at it
# par(mfrow=c(1,3))
# for (l in c(0.5:5.5)){
#   plot(df_metric$june_fcast[df_metric$lead==l], df_metric$june_obs[df_metric$lead==l], col=as.factor(df_metric$unique_model[df_metric$lead==l]))#, ylim=c(-2,2), xlim=c(-2,2))
#   abline(a=0,b=1)
#   abline(h=0,v=0)
#   plot(df_metric$jul_fcast[df_metric$lead==l], df_metric$jul_obs[df_metric$lead==l], col=as.factor(df_metric$unique_model[df_metric$lead==l]))#, ylim=c(-2,2), xlim=c(-2,2))
#   abline(a=0,b=1)
#   abline(h=0,v=0)
#   plot(df_metric$aug_fcast[df_metric$lead==l], df_metric$aug_obs[df_metric$lead==l], col=as.factor(df_metric$unique_model[df_metric$lead==l]))#, ylim=c(-2,2), xlim=c(-2,2))
#   abline(a=0,b=1)
#   abline(h=0,v=0)
# }

# #And now make plots for downscaled only
# df_metric <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/total_metric_globalforecasts.rds')
# df_metric <- df_metric[df_metric$glob_model=="downscaled",]
# df_metric$unique_model <- paste(df_metric$glob_model,df_metric$ensemble,sep="_")
# df_metric <- df_metric[df_metric$ensemble==2 |df_metric$ensemble==8 | df_metric$ensemble==10,]
# 
# ggplot(data=df_metric,aes(x=june_fcast, y=june_obs, label=year))+
#   geom_point(aes(col=unique_model),alpha=0.6)+
#   ggtitle("Downscaled: June Rule")+
#   stat_cor(method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
#   geom_hline(yintercept=0.77,col="dark grey")+
#   geom_vline(xintercept=0.77, col="dark grey")+
#   geom_text(size=2)
# 
# ggplot(data=df_metric,aes(x=jul_fcast, y=jul_obs, label=year))+
#   geom_point(aes(col=unique_model),alpha=0.6)+
#   ggtitle("June Rule")+
#   stat_cor(method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
#   geom_hline(yintercept=0.77,col="dark grey")+
#   geom_vline(xintercept=0.77, col="dark grey")+
#   geom_text(size=2)
# 
# ggplot(data=df_metric,aes(x=aug_fcast, y=aug_obs, label=year))+
#   geom_point(aes(col=unique_model),alpha=0.6)+
#   ggtitle("June Rule")+
#   stat_cor(method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
#   geom_hline(yintercept=0.77,col="dark grey")+
#   geom_vline(xintercept=0.77, col="dark grey")+
#   geom_text(size=2)
# 


#-----Step 4b: Calculate SSTA for 6-month period: apples scenario-------
#load data
obs_anom <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1b/observed_anomalies_apples.rds')
glob <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2b/global_fcasts_apples.rds')
downscaled_fcasts <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_3b/downscaled_forecasts.rds')

#prepare for merge
downscaled_fcasts$glob_model <- "downscaled"
downscaled_fcasts <- downscaled_fcasts[,c(9,1,2,3,4,5,6,7,8)]
downscaled_fcasts$ensemble <- gsub("ens","",downscaled_fcasts$ensemble)
downscaled_fcasts <- downscaled_fcasts[downscaled_fcasts$ensemble!="mean",]
downscaled_fcasts$ensemble <- as.numeric(downscaled_fcasts$ensemble)
head(downscaled_fcasts)
head(glob)
all_fcasts <- rbind(downscaled_fcasts,glob)
head(all_fcasts)

#prepare obs for merge
head(obs_anom)
obs_anom$sst_obs <- obs_anom$sst
obs_anom$sst_clim_obs <- obs_anom$sst_clim
obs_anom$sst_anom_obs <- obs_anom$obs_anom
obs_anom$sst_anom_obs_rolling <- obs_anom$roll_obs_anom
obs_anom$forecast_date <- obs_anom$date
obs_anom <- obs_anom[,c(9,10,11,12,13)]

all_data <- left_join(all_fcasts,obs_anom,"forecast_date")

globs <- c("downscaled","CanCM4i.nc","COLA-RSMAS-CCSM4.nc", "GEM-NEMO.nc",
           "GFDL-SPEAR.nc", "NASA-GEOSS2S.nc", "NCEP-CFSv2.nc")
df_metric <- as.data.frame(matrix(NA, nrow=63, ncol=10))
colnames(df_metric) <- c("glob_model","ensemble","lead", "year", "june_fcast", "june_obs", "jul_fcast", "jul_obs", "aug_fcast", "aug_obs")
counter=1
t1 <- Sys.time() #Takes 20 mins
for (g in globs){
  for (l in 0.5:5.5){ #start lead time, which then adds 6 months to get lead end time. 
    for (m in 1:24){ #1:24
      for (y in 1981:2020){  
        print(glue("model: {g}; lead {l}; ens {m}; year {y}"))
        #Get mean over Dec-may
        # get forecast dates for Jun: dec- may
        dec <- as.Date(glue('{y}-12-01'))
        may <- as.Date(glue('{y+1}-05-01')) #do it this way to make sure you get correct year
        init_date <-  ymd(glue('{y}-12-01')) %m-% months(l-0.5)
        
        jun_fcast <- mean(all_data$sst_anom[all_data$glob_model==g &
                                              all_data$ensemble==m &
                                              all_data$init_month==init_date&
                                              all_data$forecast_date>=dec & 
                                              all_data$forecast_date<=may], na.rm=T)
        jun_obs <- mean(all_data$sst_anom_obs[ all_data$glob_model==g &
                                                 all_data$ensemble==m &
                                                 all_data$init_month==init_date&
                                                 all_data$forecast_date>=dec & 
                                                 all_data$forecast_date<=may], na.rm=T)
        
        #Get mean over Jan-Jun
        init_date <-  ymd(glue('{y+1}-01-01')) %m-% months(l-0.5)
        jul_fcast <- mean(all_data$sst_anom[all_data$glob_model==g &
                                              all_data$ensemble==m &
                                              all_data$init_month==init_date&
                                              all_data$month<=6], na.rm=T)
        jul_obs <- mean(all_data$sst_anom_obs[all_data$glob_model==g &
                                                all_data$ensemble==m &
                                                all_data$init_month==init_date&
                                                all_data$month<=6], na.rm=T)
        
        #Get mean over Feb-Jul
        init_date <-  ymd(glue('{y+1}-02-01')) %m-% months(l-0.5)
        aug_fcast <- mean(all_data$sst_anom[all_data$glob_model==g &
                                              all_data$ensemble==m &
                                              all_data$init_month==init_date&
                                              all_data$month<=7 & all_data$month>=2],na.rm=T)
        aug_obs <- mean(all_data$sst_anom_obs[all_data$glob_model==g &
                                                all_data$ensemble==m &
                                                all_data$init_month==init_date&
                                                all_data$month<=7 & all_data$month>=2],na.rm=T)
        #write out data
        df_metric[counter,1] <- g
        df_metric[counter,2] <- m
        df_metric[counter,3] <- l
        df_metric[counter,4] <- y+1
        df_metric[counter,5] <- jun_fcast
        df_metric[counter,6] <- jun_obs
        df_metric[counter,7] <- jul_fcast
        df_metric[counter,8] <- jul_obs
        df_metric[counter,9] <- aug_fcast
        df_metric[counter,10] <- aug_obs
        counter=counter+1
      }
    }
  }
}
t2 <- Sys.time()
t2-t1
df_metric

saveRDS(df_metric,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_4b/total_metric_allforecasts_1980-2010_apples.rds')



ggplot(data=df_metric,aes(x=june_fcast, y=june_obs))+
  geom_point(aes(col=glob_model),alpha=0.6)+
  ggtitle("June Rule")+
  stat_cor(method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.77,col="dark grey")+
  geom_vline(xintercept=0.77, col="dark grey")+
  facet_wrap(~lead)

ggplot(data=df_metric,aes(x=jul_fcast, y=jul_obs))+
  geom_point(aes(col=glob_model),alpha=0.6)+
  ggtitle("July Rule")+
  stat_cor(method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.77,col="dark grey")+
  geom_vline(xintercept=0.77, col="dark grey")+
  facet_wrap(~lead)

ggplot(data=df_metric,aes(x=aug_fcast, y=aug_obs))+
  geom_point(aes(col=glob_model),alpha=0.6)+
  ggtitle("August Rule")+
  stat_cor(method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.77,col="dark grey")+
  geom_vline(xintercept=0.77, col="dark grey")+
  facet_wrap(~lead)

#And now make plots for downscaled only
df_metric <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/total_metric_globalforecasts.rds')
df_metric <- df_metric[df_metric$glob_model=="downscaled",]
df_metric$unique_model <- paste(df_metric$glob_model,df_metric$ensemble,sep="_")
df_metric <- df_metric[df_metric$ensemble==2 |df_metric$ensemble==8 | df_metric$ensemble==10,]

ggplot(data=df_metric,aes(x=june_fcast, y=june_obs, label=year))+
  geom_point(aes(col=unique_model),alpha=0.6)+
  ggtitle("Downscaled: June Rule")+
  stat_cor(method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.77,col="dark grey")+
  geom_vline(xintercept=0.77, col="dark grey")+
  geom_text(size=2)

ggplot(data=df_metric,aes(x=jul_fcast, y=jul_obs, label=year))+
  geom_point(aes(col=unique_model),alpha=0.6)+
  ggtitle("June Rule")+
  stat_cor(method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.77,col="dark grey")+
  geom_vline(xintercept=0.77, col="dark grey")+
  geom_text(size=2)

ggplot(data=df_metric,aes(x=aug_fcast, y=aug_obs, label=year))+
  geom_point(aes(col=unique_model),alpha=0.6)+
  ggtitle("June Rule")+
  stat_cor(method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.77,col="dark grey")+
  geom_vline(xintercept=0.77, col="dark grey")+
  geom_text(size=2)

#------Step 5a: make ensemble forecast: best case-----
#Load data
all_data <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_4a/total_metric_globalforecasts_1980-2020_bestcase.rds')
all_data <- na.omit(all_data)

#ensemble
all_data_ens <- all_data %>% group_by(lead,year) %>%  summarise_at(c("june_fcast", "june_obs",
                                                             "jul_fcast", "jul_obs",
                                                             "aug_fcast", "aug_obs"), mean, na.rm=T)
all_data_ens <- as.data.frame(all_data_ens)

#Bring in fcast quantiles
th <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2a/thresholds_globalforecasts_bestcase.rds')
th <- as.data.frame(th)

#Bring in observed quantiles
ob <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1a/observed_threshold_bestcase.rds')

all_data_ens_quants <- left_join(all_data_ens, th, "lead")
all_data_ens_quants$observed_quantile <- ob$mean

saveRDS(all_data_ens_quants,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_5a/total_metric_globalforecasts_1980-2020_bestcase_ensemble.rds')


#Plots
ggplot(data=all_data_ens_quants)+
  geom_point(aes(x=june_fcast, y=june_obs))+
  ggtitle("June Rule")+
  stat_cor(aes(x=june_fcast, y=june_obs),method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.41,col="dark grey")+
  geom_vline(aes(xintercept=ssta_quantile),col="dark grey")+
  facet_wrap(~lead)

ggplot(data=all_data_ens_quants)+
  geom_point(aes(x=jul_fcast, y=jul_obs))+
  ggtitle("July Rule")+
  stat_cor(aes(x=jul_fcast, y=jul_obs), method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.41,col="dark grey")+
  geom_vline(aes(xintercept=ssta_quantile),col="dark grey")+
  facet_wrap(~lead)

ggplot(data=all_data_ens_quants)+
  geom_point(aes(x=aug_fcast, y=aug_obs))+
  ggtitle("August Rule")+
  stat_cor(aes(x=aug_fcast, y=aug_obs),method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.41,col="dark grey")+
  geom_vline(aes(xintercept=ssta_quantile),col="dark grey")+
  facet_wrap(~lead)


#OUTDATED CODE
#I don't know what this was for????
# df_metric <- as.data.frame(matrix(NA, nrow=10, ncol=8))
# colnames(df_metric) <- c("lead", "year", "june_fcast", "june_obs", "jul_fcast", "jul_obs", "aug_fcast", "aug_obs")
# counter=1
# t1 <- Sys.time() 
# for (l in 0.5:5.5){ #start lead time, which then adds 6 months to get lead end time. 
#   for (y in 1981:2020){  
#     print(glue("lead {l}; year {y}"))
#     #Get mean over Dec-may
#     # get forecast dates for Jun: dec- may
#     dec <- as.Date(glue('{y}-12-01'))
#     may <- as.Date(glue('{y+1}-05-01')) #do it this way to make sure you get correct year
#     init_date <-  ymd(glue('{y}-12-01')) %m-% months(l-0.5)
#     
#     jun_fcast <- ifelse(length(all_data$sst_anom[
#       all_data$init_month==init_date&
#         all_data$forecast_date>=dec & 
#         all_data$forecast_date<=may])==6,
#       mean(all_data$sst_anom[
#         all_data$init_month==init_date&
#           all_data$forecast_date>=dec & 
#           all_data$forecast_date<=may], na.rm=T),
#       NA)
#     jun_obs <- mean(all_data$sst_anom_obs[
#       all_data$init_month==init_date&
#         all_data$forecast_date>=dec & 
#         all_data$forecast_date<=may], na.rm=T)
#     
#     #Get mean over Jan-Jun
#     init_date <-  ymd(glue('{y+1}-01-01')) %m-% months(l-0.5)
#     
#     jul_fcast <- ifelse(length(all_data$sst_anom[
#       all_data$init_month==init_date&
#         all_data$month<=6])==6,
#       mean(all_data$sst_anom[
#         all_data$init_month==init_date&
#           all_data$month<=6], na.rm=T),
#       NA)
#     jul_obs <- mean(all_data$sst_anom_obs_rolling[
#       all_data$init_month==init_date&
#         all_data$month<=6], na.rm=T)
#     
#     #Get mean over Feb-Jul
#     init_date <-  ymd(glue('{y+1}-02-01')) %m-% months(l-0.5)
#     aug_fcast <- ifelse(length(all_data$sst_anom[
#       all_data$init_month==init_date&
#         all_data$month<=7 & all_data$month>=2])==6,
#       mean(all_data$sst_anom[
#         all_data$init_month==init_date&
#           all_data$month<=7 & all_data$month>=2], na.rm=T),
#       NA)
#     
#     aug_obs <- mean(all_data$sst_anom_obs[
#       all_data$init_month==init_date&
#         all_data$month<=7 & all_data$month>=2],na.rm=T)
#     #write out data
#     df_metric[counter,1] <- l
#     df_metric[counter,2] <- y+1
#     df_metric[counter,3] <- jun_fcast
#     df_metric[counter,4] <- jun_obs
#     df_metric[counter,5] <- jul_fcast
#     df_metric[counter,6] <- jul_obs
#     df_metric[counter,7] <- aug_fcast
#     df_metric[counter,8] <- aug_obs
#     counter=counter+1
#   }
# }
# t2 <- Sys.time()
# t2-t1
# df_metric
# 
# saveRDS(df_metric,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/total_metric_globalforecasts_allyears_ensemble.rds')

# 
# ##PLOTS
# df_metric <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/total_metric_globalforecasts_allyears_ensemble.rds')
# #bring in quantiles
# quants <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/ssta_total_forecast_quantiles.rds')
# #merge
# df_metric <- left_join(df_metric, quants, "lead")
# 




#------Step 5b: make ensemble forecast: apples-to-apples-----
#Load data
all_data <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_4b/total_metric_allforecasts_1980-2010_apples.rds')
all_data <- na.omit(all_data)

#ensemble
all_data_ens <- all_data %>% group_by(lead,year) %>%  summarise_at(c("june_fcast", "june_obs",
                                                                     "jul_fcast", "jul_obs",
                                                                     "aug_fcast", "aug_obs"), mean, na.rm=T)

#Bring in previous quantiles
th <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2b/thresholds_globalforecasts_apples.rds')

#Bring in observed quantiles
ob <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1b/observed_threshold_apples.rds')

all_data_ens_quants <- left_join(all_data_ens, th, "lead")
all_data_ens_quants$observed_quantile <- ob


saveRDS(all_data_ens_quants,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_5b/total_metric_globalforecasts_1980-2010_apples_ensemble.rds')


#------Step 5c: make ensemble forecast: only-3-----
#Load data
all_data <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_4b/total_metric_allforecasts_1980-2010_apples.rds')
all_data <- na.omit(all_data)

all_data <- all_data[all_data$glob_model=="CanCM4i.nc",]
all_data <- all_data[all_data$ensemble==2 | all_data$ensemble==8 | all_data$ensemble==10,]

#ensemble
all_data_ens <- all_data %>% group_by(lead,year) %>%  summarise_at(c("june_fcast", "june_obs",
                                                                     "jul_fcast", "jul_obs",
                                                                     "aug_fcast", "aug_obs"), mean, na.rm=T)

#Bring in previous quantiles
th <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2c/thresholds_globalforecasts_only3.rds')

#Bring in observed quantiles
ob <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1b/observed_threshold_apples.rds')

all_data_ens_quants <- left_join(all_data_ens, th, "lead")
all_data_ens_quants$observed_quantile <- ob


saveRDS(all_data_ens_quants,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_5c/total_metric_globalforecasts_1980-2010_only3_ensemble.rds')


#------Step 6b: Calculate SSTA for 6-month period: downscaled forecasts apples-----
# Jun = Dec-May mean
# Jul = Jan-Jun mean
# Aug = Feb=July mean

all_data <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_3b/downscaled_forecasts.rds')
obs <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1b/observed_anomalies_apples.rds')
obs$forecast_date <- obs$date
all_data <- left_join(all_data,obs[,8:9],"forecast_date")
all_data$year <- lubridate::year(all_data$forecast_date)
# all_data$month <- lubridate::month(all_data$forecast_date)
# all_data <- all_data[all_data$glob_model!="downscaled",]
all_data <- na.omit(all_data)

models <- c("ens10", "ens2", "ens8")
df_metric <- as.data.frame(matrix(NA, nrow=1, ncol=9))
colnames(df_metric) <- c("model","lead", "year", "june_fcast", "june_obs", "jul_fcast", "jul_obs", "aug_fcast", "aug_obs")
counter=1
t1 <- Sys.time() #Takes 12 mins
for (m in models){
  for (l in 0.5:5.5){ #start lead time, which then adds 6 months to get lead end time.
    for (y in 1982:2010){
      print(glue("lead {l}; year {y}"))
      #Get mean over Dec-may
      # get forecast dates for Jun: dec- may
      dec <- as.Date(glue('{y}-12-01'))
      may <- as.Date(glue('{y+1}-05-01')) #do it this way to make sure you get correct year
      init_date <-  ymd(glue('{y}-12-01')) %m-% months(l-0.5)
      
      jun_fcast <- ifelse(length(all_data$sst_anom[all_data$ensemble==m & all_data$init_month==init_date &
                                                     all_data$forecast_date>=dec & all_data$forecast_date<=may])==6,
                          mean(all_data$sst_anom[all_data$ensemble==m & all_data$init_month==init_date & all_data$forecast_date>=dec &
                                                   all_data$forecast_date<=may], na.rm=T), NA)
      
      
      jun_obs <- mean(all_data$roll_obs_anom[ all_data$init_month==init_date&
                                              all_data$forecast_date>=dec &
                                              all_data$forecast_date<=may], na.rm=T)
      
      #Get mean over Jan-Jun
      init_date <-  ymd(glue('{y+1}-01-01')) %m-% months(l-0.5)
      jul_fcast <- ifelse(length(all_data$sst_anom[all_data$ensemble==m & all_data$init_month==init_date &
                                                     all_data$month<=6])==6,
                          mean(all_data$sst_anom[all_data$ensemble==m & all_data$init_month==init_date &
                                                   all_data$month<=6], na.rm=T), NA)
      
      jul_obs <- mean(all_data$roll_obs_anom[ all_data$init_month==init_date&
                                                      all_data$month<=6], na.rm=T)
      
      #Get mean over Feb-Jul
      init_date <-  ymd(glue('{y+1}-02-01')) %m-% months(l-0.5)
      aug_fcast <- ifelse(length(all_data$sst_anom[all_data$ensemble==m & all_data$init_month==init_date &
                                                     all_data$month<=7 & all_data$month>=2])==6,
                          mean(all_data$sst_anom[all_data$ensemble==m & all_data$init_month==init_date &
                                                   all_data$month<=7 & all_data$month>=2], na.rm=T), NA)
      
      aug_obs <- mean(all_data$roll_obs_anom[ all_data$init_month==init_date&
                                              all_data$month<=7 & all_data$month>=2],na.rm=T)
      #write out data
      df_metric[counter,1] <- m
      df_metric[counter,2] <- l
      df_metric[counter,3] <- y+1
      df_metric[counter,4] <- jun_fcast
      df_metric[counter,5] <- jun_obs
      df_metric[counter,6] <- jul_fcast
      df_metric[counter,7] <- jul_obs
      df_metric[counter,8] <- aug_fcast
      df_metric[counter,9] <- aug_obs
      counter=counter+1
    }
  }
}
t2 <- Sys.time()
t2-t1
df_metric
saveRDS(df_metric,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_6b/total_metric_downscaledforecasts_ens_1982-2010.rds')


#PLOTS
ggplot(data=df_metric)+
  geom_point(aes(x=june_fcast, y=june_obs))+
  ggtitle("June Rule")+
  stat_cor(aes(x=june_fcast, y=june_obs),method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.49,col="dark grey")+
  # geom_vline(aes(xintercept=ssta_quantile), col="dark grey")+
  facet_wrap(~lead)

ggplot(data=df_metric)+
  geom_point(aes(x=jul_fcast, y=jul_obs))+
  ggtitle("July Rule")+
  stat_cor(aes(x=jul_fcast, y=jul_obs),method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.49,col="dark grey")+
  # geom_vline(aes(xintercept=ssta_quantile), col="dark grey")+
  facet_wrap(~lead)

ggplot(data=df_metric)+
  geom_point(aes(x=aug_fcast, y=aug_obs))+
  ggtitle("August Rule")+
  stat_cor(aes(x=aug_fcast, y=aug_obs),method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.49,col="dark grey")+
  # geom_vline(aes(xintercept=ssta_quantile), col="dark grey")+
  facet_wrap(~lead)



#------Step 7b: make ensemble downscaled forecast: apples------
#Load data
all_data <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_6b/total_metric_downscaledforecasts_ens_1982-2010.rds')

#ensemble
all_data_ens <- all_data %>% group_by(lead,year) %>%  summarise_at(c("june_fcast", "june_obs",
                                                                     "jul_fcast", "jul_obs",
                                                                     "aug_fcast", "aug_obs"), mean, na.rm=T)
all_data_ens <- as.data.frame(all_data_ens)

#Bring in previous quantiles
th <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_3b/ssta_total_forecast_quantiles_downscaled.rds')
th <- as.data.frame(th)

#Bring in observed quantiles
ob <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1b/observed_threshold_apples.rds')

all_data_ens_quants <- left_join(all_data_ens, th, "lead")
all_data_ens_quants$observed_quantile <- ob

saveRDS(all_data_ens_quants,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_7b/total_metric_globalforecasts_1980-2010_apples_ensemble.rds')




#Plots
g1 <- ggplot(data=all_data_ens_quants)+
  geom_point(aes(x=june_fcast, y=june_obs))+
  ggtitle("June Rule")+
  stat_cor(aes(x=june_fcast, y=june_obs),method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.64,col="dark grey")+
  geom_vline(xintercept=0.36,col="dark grey")+
  lims(x=c(-1.55,1.35),y=c(-1.5,1.7))

g2 <- ggplot(data=all_data_ens_quants)+
  geom_point(aes(x=jul_fcast, y=jul_obs))+
  ggtitle("July Rule")+
  stat_cor(aes(x=jul_fcast, y=jul_obs), method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.64,col="dark grey")+
  geom_vline(xintercept=0.45,col="dark grey")+
  lims(x=c(-1.55,1.35),y=c(-1.5,1.7))

g3 <- ggplot(data=all_data_ens_quants)+
  geom_point(aes(x=aug_fcast, y=aug_obs))+
  ggtitle("August Rule")+
  stat_cor(aes(x=aug_fcast, y=aug_obs),method = "pearson", label.x = -1.5, label.y = 1.8, size=3)+
  geom_hline(yintercept=0.64,col="dark grey")+
  geom_vline(xintercept=0.44,col="dark grey")+
  lims(x=c(-1.55,1.35),y=c(-1.5,1.7))

grid.arrange(g1, g2, g3, nrow=1)


#-----Step 8a: Estimate Skill: best case scenario------
#TNR, TPR, ACC, SEDI, and Cor

df_metric <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_5a/total_metric_globalforecasts_1980-2020_bestcase_ensemble.rds')

metric_summary <- as.data.frame(matrix(NA,nrow=4,ncol=22))
colnames(metric_summary) <- c("lead",
                              "cor_jun","cor_jul","cor_aug",
                              "cor_sig_jun","cor_sig_jul","cor_sig_aug",
                              "tpr_jun","tpr_jul","tpr_aug",
                              "tnr_jun","tnr_jul","tnr_aug",
                              "acc_jun","acc_jul","acc_aug",
                              "sedi_jun","sedi_jul","sedi_aug",
                              "bss_jun","bss_jul","bss_aug")
q = unique(df_metric$observed_quantile)
counter=1
for (l in 0.5:5.5){ #start lead time, which then adds 6 months to get lead end time. 
  
  #Correlation coefficients
  # cor_jun <- cor(df_metric$june_fcast[df_metric$lead==l], df_metric$june_obs[df_metric$lead==l], use="na.or.complete")
  # cor_jul <- cor(df_metric$jul_fcast[df_metric$lead==l], df_metric$jul_obs[ df_metric$lead==l], use="na.or.complete")
  # cor_aug <- cor(df_metric$aug_fcast[ df_metric$lead==l], df_metric$aug_obs[ df_metric$lead==l], use="na.or.complete")
  
  data_temp <- df_metric[df_metric$lead==l,]
  
  #Correlation coefficients
  cor_i_jun <- cor(data_temp$june_fcast, data_temp$june_obs, use="na.or.complete")
  cor_i_jul <- cor(data_temp$jul_fcast, data_temp$jul_obs, use="na.or.complete")
  cor_i_aug <- cor(data_temp$aug_fcast, data_temp$aug_obs, use="na.or.complete")
  
  #Calculate N effective degrees of freedom to determine if correlation is significant 
  N <- length(unique(data_temp$year))
  tau <- 0:(N-1)
  dm <- as.data.frame(as.matrix(NA,nrow=N,ncol=3))
  for (t in 1:length(tau)){
    if(t<N-1){
      rx_jun <- cor(data_temp$june_fcast[1:(N-tau[t])],data_temp$june_fcast[(tau[t]+1):N])
      ry_jun <- cor(data_temp$june_obs[1:(N-tau[t])],data_temp$june_obs[(tau[t]+1):N], use="na.or.complete")
      
      rx_jul <- cor(data_temp$jul_fcast[1:(N-tau[t])],data_temp$jul_fcast[(tau[t]+1):N])
      ry_jul <- cor(data_temp$jul_obs[1:(N-tau[t])],data_temp$jul_obs[(tau[t]+1):N], use="na.or.complete")
      
      rx_aug <- cor(data_temp$aug_fcast[1:(N-tau[t])],data_temp$aug_fcast[(tau[t]+1):N])
      ry_aug <- cor(data_temp$aug_obs[1:(N-tau[t])],data_temp$aug_obs[(tau[t]+1):N], use="na.or.complete")
    } else {
      rx_jun <- rx_jul <- rx_aug <- 1
      ry_jun <- ry_jul <- ry_aug <- 1
    }
    dm[t,1] <- (1-tau[t]/N) * rx_jun * ry_jun
    dm[t,2] <- (1-tau[t]/N) * rx_jul * ry_jul
    dm[t,3] <- (1-tau[t]/N) * rx_aug * ry_aug
  }
  NEFF_jun <- N/sum(dm$V1, na.rm=T)  
  NEFF_jul <- N/sum(dm$V2, na.rm=T)  
  NEFF_aug <- N/sum(dm$V3, na.rm=T)  
  ci_jun <- CorCI(cor_i_jun,n=NEFF_jun, conf.level = 0.95)
  ci_jul <- CorCI(cor_i_jul,n=NEFF_jul, conf.level = 0.95)
  ci_aug <- CorCI(cor_i_aug,n=NEFF_aug, conf.level = 0.95)
  ci_l_jun <- ifelse(ci_jun[2]>0,1,0 )
  ci_l_jul <- ifelse(ci_jul[2]>0,1,0 )
  ci_l_aug <- ifelse(ci_aug[2]>0,1,0 )
  
  #TPR: TP/TP+FN (not x has to be observed)
  tpr_jun <- tpr(o = df_metric$june_obs[ df_metric$lead==l], f = df_metric$june_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tpr_jul <- tpr(o = df_metric$jul_obs[ df_metric$lead==l], f = df_metric$jul_fcast[df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tpr_aug <- tpr(o = df_metric$aug_obs[df_metric$lead==l], f = df_metric$aug_fcast[df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  
  #TNR: TN/TN+FP
  tnr_jun <- tnr(o = df_metric$june_obs[ df_metric$lead==l], f = df_metric$june_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tnr_jul <- tnr(o = df_metric$jul_obs[df_metric$lead==l], f = df_metric$jul_fcast[df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tnr_aug <- tnr(o = df_metric$aug_obs[ df_metric$lead==l], f = df_metric$aug_fcast [df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  
  #ACC: TN + TN / total
  acc_jun <- acc(o = df_metric$june_obs[ df_metric$lead==l], f = df_metric$june_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  acc_jul <- acc(o = df_metric$jul_obs[df_metric$lead==l], f = df_metric$jul_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  acc_aug <- acc(o = df_metric$aug_obs[ df_metric$lead==l], f = df_metric$aug_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  
  obs_jun <- df_metric$june_obs[ df_metric$lead==l]
  pred_jun <- df_metric$june_fcast[ df_metric$lead==l]
  obs_jun <- obs_jun[!is.na(obs_jun)]
  pred_jun <- pred_jun[!is.na(pred_jun)]
  v_jun <- verify(obs = ifelse(obs_jun>q,1,0), pred = ifelse(pred_jun>unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0), frcst.type = "binary", obs.type = "binary")
  
  obs_jul <- df_metric$jul_obs[ df_metric$lead==l]
  pred_jul <- df_metric$jul_fcast[df_metric$lead==l]
  obs_jul <- obs_jul[!is.na(obs_jul)]
  pred_jul <- pred_jul[!is.na(pred_jul)]
  v_jul <- verify(obs = ifelse(obs_jul>q,1,0), pred = ifelse(pred_jul>unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0), frcst.type = "binary", obs.type = "binary")
  
  obs_aug <- df_metric$aug_obs[df_metric$lead==l]
  pred_aug <- df_metric$aug_fcast[ df_metric$lead==l]
  obs_aug <- obs_aug[!is.na(obs_aug)]
  pred_aug <- pred_aug[!is.na(pred_aug)]
  v_aug <- verify(obs = ifelse(obs_aug>q,1,0), pred = ifelse(pred_aug>unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0), frcst.type = "binary", obs.type = "binary")
  
  #BS and BSS
  obs_jun <- ifelse(obs_jun <= q,1,0)
  pred_jun <- ifelse(pred_jun <= unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0)
  BS_jun <- mean((pred_jun - obs_jun)^2)
  BS_ref_jun <-  0.74
  BSS_jun <- 1 - (BS_jun/BS_ref_jun)
  
  obs_jul <- ifelse(obs_jul <= q,1,0)
  pred_jul <- ifelse(pred_jul <= unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0)
  BS_jul <- mean((pred_jul - obs_jul)^2)
  BS_ref_jul <-  0.74
  BSS_jul <- 1 - (BS_jul/BS_ref_jul)
  
  obs_aug <- ifelse(obs_aug <= q,1,0)
  pred_aug <- ifelse(pred_aug <= unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0)
  BS_aug <- mean((pred_aug - obs_aug)^2)
  BS_ref_aug <-  0.74
  BSS_aug <- 1 - (BS_aug/BS_ref_aug)
  
  metric_summary[counter,1] <- l
  metric_summary[counter,2] <- cor_i_jun
  metric_summary[counter,3] <- cor_i_jul
  metric_summary[counter,4] <- cor_i_aug
  metric_summary[counter,5] <- ci_l_jun
  metric_summary[counter,6] <- ci_l_jul
  metric_summary[counter,7] <- ci_l_aug
  metric_summary[counter,8] <- tpr_jun
  metric_summary[counter,9] <- tpr_jul
  metric_summary[counter,10] <- tpr_aug
  metric_summary[counter,11] <- tnr_jun
  metric_summary[counter,12] <- tnr_jul
  metric_summary[counter,13] <- tnr_aug
  metric_summary[counter,14] <- acc_jun
  metric_summary[counter,15] <- acc_jul
  metric_summary[counter,16] <- acc_aug
  metric_summary[counter,17] <- v_jun$SEDI
  metric_summary[counter,18] <- v_jul$SEDI
  metric_summary[counter,19] <- v_aug$SEDI
  metric_summary[counter,20] <- BSS_jun
  metric_summary[counter,21] <- BSS_jul
  metric_summary[counter,22] <- BSS_aug
  counter=counter+1
  
}

metric_summary2 <- metric_summary %>% pivot_longer(c(2:4), names_to = "correlation_months", values_to="correlation_values") %>% 
  pivot_longer(c(2:4), names_to = "sig_months", values_to="sig_values") %>% 
  pivot_longer(c(2:4), names_to = "tpr_months", values_to="tpr_values") %>% 
  pivot_longer(c(2:4), names_to = "tnr_months", values_to="tnr_values") %>% 
  pivot_longer(c(2:4), names_to = "acc_months", values_to="acc_values") %>% 
  pivot_longer(c(2:4), names_to = "sedi_months", values_to="sedi_values") %>% 
  pivot_longer(c(2:4), names_to = "bss_months", values_to="bss_values")
metric_summary2 <- metric_summary2 %>%  separate(correlation_months, c(NA, "correlation_months"))
metric_summary2 <- metric_summary2 %>%  separate(sig_months, c(NA,NA, "cor_sig_months"))
metric_summary2 <- metric_summary2 %>%  separate(tpr_months, c(NA, "tpr_months"))
metric_summary2 <- metric_summary2 %>%  separate(tnr_months, c(NA, "tnr_months"))
metric_summary2 <- metric_summary2 %>%  separate(acc_months, c(NA, "acc_months"))
metric_summary2 <- metric_summary2 %>%  separate(sedi_months, c(NA, "sedi_months"))
metric_summary2 <- metric_summary2 %>%  separate(bss_months, c(NA, "bss_months"))

saveRDS(metric_summary2,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_8a/metric_summary_global_1980-2020_bestcasescenario.rds')

#Correlation plots
metric_summary2$correlation_months <- as.factor(metric_summary2$correlation_months)
metric_summary2$correlation_months <- factor(metric_summary2$correlation_months, levels=rev(levels(metric_summary2$correlation_months)))
cor_plot <- ggplot(data = metric_summary2, aes(x=correlation_months,y=lead,fill=correlation_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(correlation_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.31, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  ggtitle("Correlation")

#TPR plots
metric_summary2$tpr_months <- as.factor(metric_summary2$tpr_months)
metric_summary2$tpr_months <- factor(metric_summary2$tpr_months, levels=rev(levels(metric_summary2$tpr_months)))
tpr_plot <- ggplot(data = metric_summary2, aes(x=tpr_months,y=lead,fill=tpr_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tpr_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TPR") +
  ggtitle("True Positive Rate")

#TNR
metric_summary2$tnr_months <- as.factor(metric_summary2$tnr_months)
metric_summary2$tnr_months <- factor(metric_summary2$tnr_months, levels=rev(levels(metric_summary2$tnr_months)))
tnr_plot <- ggplot(data = metric_summary2, aes(x=tnr_months,y=lead,fill=tnr_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tnr_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TNR") +
  ggtitle("True Negative Rate")

#ACC
metric_summary2$acc_months <- as.factor(metric_summary2$acc_months)
metric_summary2$acc_months <- factor(metric_summary2$acc_months, levels=rev(levels(metric_summary2$acc_months)))
acc_plot <- ggplot(data = metric_summary2, aes(x=acc_months,y=lead,fill=acc_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(acc_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.62, limit = c(0,1), space = "Lab",
                       name="Accuracy") +
  ggtitle("Accuracy")
# (0.26*0.26) + (0.74*0.74) = 0.6152

#SEDI
metric_summary2$sedi_months <- as.factor(metric_summary2$sedi_months)
metric_summary2$sedi_months <- factor(metric_summary2$sedi_months, levels=rev(levels(metric_summary2$sedi_months)))
sedi_plot <- ggplot(data = metric_summary2, aes(x=sedi_months,y=lead,fill=sedi_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(sedi_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="SEDI") +
  ggtitle("SEDI")

#BSS
metric_summary2$bss_months <- as.factor(metric_summary2$bss_months)
metric_summary2$bss_months <- factor(metric_summary2$bss_months, levels=rev(levels(metric_summary2$bss_months)))
bss_plot <- ggplot(data = metric_summary2, aes(x=bss_months,y=lead,fill=bss_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(bss_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="BSS") +
  ggtitle("Brier Skill Score")

tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/fig_skill_plots_global_bestcase_step8a.tiff', res=200, units="in",width=12,height=8)
grid.arrange(cor_plot, tpr_plot, tnr_plot, acc_plot, sedi_plot, bss_plot, nrow=2)
dev.off()


#-----Step 8b: Estimate Skill: apples scenario------
#TNR, TPR, ACC, SEDI, and Cor

df_metric <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_5b/total_metric_globalforecasts_1980-2010_apples_ensemble.rds')

metric_summary <- as.data.frame(matrix(NA,nrow=4,ncol=19))
colnames(metric_summary) <- c("lead",
                              "cor_jun","cor_jul","cor_aug",
                              "tpr_jun","tpr_jul","tpr_aug",
                              "tnr_jun","tnr_jul","tnr_aug",
                              "acc_jun","acc_jul","acc_aug",
                              "sedi_jun","sedi_jul","sedi_aug",
                              "bss_jun","bss_jul","bss_aug")
q = unique(df_metric$observed_quantile)
counter=1
for (l in 0.5:5.5){ #start lead time, which then adds 6 months to get lead end time. 
  
  #Correlation coefficients
  cor_jun <- cor(df_metric$june_fcast[df_metric$lead==l], df_metric$june_obs[df_metric$lead==l], use="na.or.complete")
  cor_jul <- cor(df_metric$jul_fcast[df_metric$lead==l], df_metric$jul_obs[ df_metric$lead==l], use="na.or.complete")
  cor_aug <- cor(df_metric$aug_fcast[ df_metric$lead==l], df_metric$aug_obs[ df_metric$lead==l], use="na.or.complete")
  
  #TPR: TP/TP+FN (not x has to be observed)
  tpr_jun <- tpr(o = df_metric$june_obs[ df_metric$lead==l], f = df_metric$june_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tpr_jul <- tpr(o = df_metric$jul_obs[ df_metric$lead==l], f = df_metric$jul_fcast[df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tpr_aug <- tpr(o = df_metric$aug_obs[df_metric$lead==l], f = df_metric$aug_fcast[df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  
  #TNR: TN/TN+FP
  tnr_jun <- tnr(o = df_metric$june_obs[ df_metric$lead==l], f = df_metric$june_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tnr_jul <- tnr(o = df_metric$jul_obs[df_metric$lead==l], f = df_metric$jul_fcast[df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tnr_aug <- tnr(o = df_metric$aug_obs[ df_metric$lead==l], f = df_metric$aug_fcast [df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  
  #ACC: TN + TN / total
  acc_jun <- acc(o = df_metric$june_obs[ df_metric$lead==l], f = df_metric$june_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  acc_jul <- acc(o = df_metric$jul_obs[df_metric$lead==l], f = df_metric$jul_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  acc_aug <- acc(o = df_metric$aug_obs[ df_metric$lead==l], f = df_metric$aug_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  
  obs_jun <- df_metric$june_obs[ df_metric$lead==l]
  pred_jun <- df_metric$june_fcast[ df_metric$lead==l]
  obs_jun <- obs_jun[!is.na(obs_jun)]
  pred_jun <- pred_jun[!is.na(pred_jun)]
  v_jun <- verify(obs = ifelse(obs_jun>q,1,0), pred = ifelse(pred_jun>unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0), frcst.type = "binary", obs.type = "binary")
  
  obs_jul <- df_metric$jul_obs[ df_metric$lead==l]
  pred_jul <- df_metric$jul_fcast[df_metric$lead==l]
  obs_jul <- obs_jul[!is.na(obs_jul)]
  pred_jul <- pred_jul[!is.na(pred_jul)]
  v_jul <- verify(obs = ifelse(obs_jul>q,1,0), pred = ifelse(pred_jul>unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0), frcst.type = "binary", obs.type = "binary")
  
  obs_aug <- df_metric$aug_obs[df_metric$lead==l]
  pred_aug <- df_metric$aug_fcast[ df_metric$lead==l]
  obs_aug <- obs_aug[!is.na(obs_aug)]
  pred_aug <- pred_aug[!is.na(pred_aug)]
  v_aug <- verify(obs = ifelse(obs_aug>q,1,0), pred = ifelse(pred_aug>unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0), frcst.type = "binary", obs.type = "binary")
  
  #BS and BSS
  obs_jun <- ifelse(obs_jun <= q,1,0)
  pred_jun <- ifelse(pred_jun <= unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0)
  BS_jun <- mean((pred_jun - obs_jun)^2)
  BS_ref_jun <-  0.82
  BSS_jun <- 1 - (BS_jun/BS_ref_jun)
  
  obs_jul <- ifelse(obs_jul <= q,1,0)
  pred_jul <- ifelse(pred_jul <= unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0)
  BS_jul <- mean((pred_jul - obs_jul)^2)
  BS_ref_jul <-  0.82
  BSS_jul <- 1 - (BS_jul/BS_ref_jul)
  
  obs_aug <- ifelse(obs_aug <= q,1,0)
  pred_aug <- ifelse(pred_aug <= unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0)
  BS_aug <- mean((pred_aug - obs_aug)^2)
  BS_ref_aug <-  0.82
  BSS_aug <- 1 - (BS_aug/BS_ref_aug)
  
  metric_summary[counter,1] <- l
  metric_summary[counter,2] <- cor_jun
  metric_summary[counter,3] <- cor_jul
  metric_summary[counter,4] <- cor_aug
  metric_summary[counter,5] <- tpr_jun
  metric_summary[counter,6] <- tpr_jul
  metric_summary[counter,7] <- tpr_aug
  metric_summary[counter,8] <- tnr_jun
  metric_summary[counter,9] <- tnr_jul
  metric_summary[counter,10] <- tnr_aug
  metric_summary[counter,11] <- acc_jun
  metric_summary[counter,12] <- acc_jul
  metric_summary[counter,13] <- acc_aug
  metric_summary[counter,14] <- v_jun$SEDI
  metric_summary[counter,15] <- v_jul$SEDI
  metric_summary[counter,16] <- v_aug$SEDI
  metric_summary[counter,17] <- BSS_jun
  metric_summary[counter,18] <- BSS_jul
  metric_summary[counter,19] <- BSS_aug
  counter=counter+1
  
}

metric_summary2 <- metric_summary %>% pivot_longer(c(2:4), names_to = "correlation_months", values_to="correlation_values") %>% 
  pivot_longer(c(2:4), names_to = "tpr_months", values_to="tpr_values") %>% 
  pivot_longer(c(2:4), names_to = "tnr_months", values_to="tnr_values") %>% 
  pivot_longer(c(2:4), names_to = "acc_months", values_to="acc_values") %>% 
  pivot_longer(c(2:4), names_to = "sedi_months", values_to="sedi_values") %>% 
  pivot_longer(c(2:4), names_to = "bss_months", values_to="bss_values")
metric_summary2 <- metric_summary2 %>%  separate(correlation_months, c(NA, "correlation_months"))
metric_summary2 <- metric_summary2 %>%  separate(tpr_months, c(NA, "tpr_months"))
metric_summary2 <- metric_summary2 %>%  separate(tnr_months, c(NA, "tnr_months"))
metric_summary2 <- metric_summary2 %>%  separate(acc_months, c(NA, "acc_months"))
metric_summary2 <- metric_summary2 %>%  separate(sedi_months, c(NA, "sedi_months"))
metric_summary2 <- metric_summary2 %>%  separate(bss_months, c(NA, "bss_months"))

saveRDS(metric_summary2,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_8b/metric_summary_global_1980-2010_applesscenario.rds')

#Correlation plots
metric_summary2$correlation_months <- as.factor(metric_summary2$correlation_months)
metric_summary2$correlation_months <- factor(metric_summary2$correlation_months, levels=rev(levels(metric_summary2$correlation_months)))
cor_plot <- ggplot(data = metric_summary2, aes(x=correlation_months,y=lead,fill=correlation_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(correlation_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.36, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  ggtitle("Correlation")

#TPR plots
metric_summary2$tpr_months <- as.factor(metric_summary2$tpr_months)
metric_summary2$tpr_months <- factor(metric_summary2$tpr_months, levels=rev(levels(metric_summary2$tpr_months)))
tpr_plot <- ggplot(data = metric_summary2, aes(x=tpr_months,y=lead,fill=tpr_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tpr_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TPR") +
  ggtitle("True Positive Rate")

#TNR
metric_summary2$tnr_months <- as.factor(metric_summary2$tnr_months)
metric_summary2$tnr_months <- factor(metric_summary2$tnr_months, levels=rev(levels(metric_summary2$tnr_months)))
tnr_plot <- ggplot(data = metric_summary2, aes(x=tnr_months,y=lead,fill=tnr_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tnr_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TNR") +
  ggtitle("True Negative Rate")

#ACC
metric_summary2$acc_months <- as.factor(metric_summary2$acc_months)
metric_summary2$acc_months <- factor(metric_summary2$acc_months, levels=rev(levels(metric_summary2$acc_months)))
acc_plot <- ggplot(data = metric_summary2, aes(x=acc_months,y=lead,fill=acc_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(acc_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.71, limit = c(0,1), space = "Lab",
                       name="Accuracy") +
  ggtitle("Accuracy")
# (0.18*0.18) + (0.82*0.82) = 0.7048

#SEDI
metric_summary2$sedi_months <- as.factor(metric_summary2$sedi_months)
metric_summary2$sedi_months <- factor(metric_summary2$sedi_months, levels=rev(levels(metric_summary2$sedi_months)))
sedi_plot <- ggplot(data = metric_summary2, aes(x=sedi_months,y=lead,fill=sedi_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(sedi_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="SEDI") +
  ggtitle("SEDI")

#BSS
metric_summary2$bss_months <- as.factor(metric_summary2$bss_months)
metric_summary2$bss_months <- factor(metric_summary2$bss_months, levels=rev(levels(metric_summary2$bss_months)))
bss_plot <- ggplot(data = metric_summary2, aes(x=bss_months,y=lead,fill=bss_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(bss_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="BSS") +
  ggtitle("Brier Skill Score")

tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/fig_skill_plots_global_apples_step8b.tiff', res=200, units="in",width=12,height=8)
grid.arrange(cor_plot, tpr_plot, tnr_plot, acc_plot, sedi_plot, bss_plot, nrow=2)
dev.off()

#-----Step 8c: Estimate Skill: only 3 scenario------
#TNR, TPR, ACC, SEDI, and Cor

df_metric <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_5c/total_metric_globalforecasts_1980-2010_only3_ensemble.rds')

metric_summary <- as.data.frame(matrix(NA,nrow=4,ncol=19))
colnames(metric_summary) <- c("lead",
                              "cor_jun","cor_jul","cor_aug",
                              "tpr_jun","tpr_jul","tpr_aug",
                              "tnr_jun","tnr_jul","tnr_aug",
                              "acc_jun","acc_jul","acc_aug",
                              "sedi_jun","sedi_jul","sedi_aug",
                              "bss_jun","bss_jul","bss_aug")
q = unique(df_metric$observed_quantile)
counter=1
for (l in 0.5:5.5){ #start lead time, which then adds 6 months to get lead end time. 
  
  #Correlation coefficients
  cor_jun <- cor(df_metric$june_fcast[df_metric$lead==l], df_metric$june_obs[df_metric$lead==l], use="na.or.complete")
  cor_jul <- cor(df_metric$jul_fcast[df_metric$lead==l], df_metric$jul_obs[ df_metric$lead==l], use="na.or.complete")
  cor_aug <- cor(df_metric$aug_fcast[ df_metric$lead==l], df_metric$aug_obs[ df_metric$lead==l], use="na.or.complete")
  
  #TPR: TP/TP+FN (not x has to be observed)
  tpr_jun <- tpr(o = df_metric$june_obs[ df_metric$lead==l], f = df_metric$june_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tpr_jul <- tpr(o = df_metric$jul_obs[ df_metric$lead==l], f = df_metric$jul_fcast[df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tpr_aug <- tpr(o = df_metric$aug_obs[df_metric$lead==l], f = df_metric$aug_fcast[df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  
  #TNR: TN/TN+FP
  tnr_jun <- tnr(o = df_metric$june_obs[ df_metric$lead==l], f = df_metric$june_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tnr_jul <- tnr(o = df_metric$jul_obs[df_metric$lead==l], f = df_metric$jul_fcast[df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tnr_aug <- tnr(o = df_metric$aug_obs[ df_metric$lead==l], f = df_metric$aug_fcast [df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  
  #ACC: TN + TN / total
  acc_jun <- acc(o = df_metric$june_obs[ df_metric$lead==l], f = df_metric$june_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  acc_jul <- acc(o = df_metric$jul_obs[df_metric$lead==l], f = df_metric$jul_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  acc_aug <- acc(o = df_metric$aug_obs[ df_metric$lead==l], f = df_metric$aug_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  
  obs_jun <- df_metric$june_obs[ df_metric$lead==l]
  pred_jun <- df_metric$june_fcast[ df_metric$lead==l]
  obs_jun <- obs_jun[!is.na(obs_jun)]
  pred_jun <- pred_jun[!is.na(pred_jun)]
  v_jun <- verify(obs = ifelse(obs_jun>q,1,0), pred = ifelse(pred_jun>unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0), frcst.type = "binary", obs.type = "binary")
  
  obs_jul <- df_metric$jul_obs[ df_metric$lead==l]
  pred_jul <- df_metric$jul_fcast[df_metric$lead==l]
  obs_jul <- obs_jul[!is.na(obs_jul)]
  pred_jul <- pred_jul[!is.na(pred_jul)]
  v_jul <- verify(obs = ifelse(obs_jul>q,1,0), pred = ifelse(pred_jul>unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0), frcst.type = "binary", obs.type = "binary")
  
  obs_aug <- df_metric$aug_obs[df_metric$lead==l]
  pred_aug <- df_metric$aug_fcast[ df_metric$lead==l]
  obs_aug <- obs_aug[!is.na(obs_aug)]
  pred_aug <- pred_aug[!is.na(pred_aug)]
  v_aug <- verify(obs = ifelse(obs_aug>q,1,0), pred = ifelse(pred_aug>unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0), frcst.type = "binary", obs.type = "binary")
  
  #BS and BSS
  obs_jun <- ifelse(obs_jun <= q,1,0)
  pred_jun <- ifelse(pred_jun <= unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0)
  BS_jun <- mean((pred_jun - obs_jun)^2)
  BS_ref_jun <-  0.82
  BSS_jun <- 1 - (BS_jun/BS_ref_jun)
  
  obs_jul <- ifelse(obs_jul <= q,1,0)
  pred_jul <- ifelse(pred_jul <= unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0)
  BS_jul <- mean((pred_jul - obs_jul)^2)
  BS_ref_jul <-  0.82
  BSS_jul <- 1 - (BS_jul/BS_ref_jul)
  
  obs_aug <- ifelse(obs_aug <= q,1,0)
  pred_aug <- ifelse(pred_aug <= unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0)
  BS_aug <- mean((pred_aug - obs_aug)^2)
  BS_ref_aug <-  0.82
  BSS_aug <- 1 - (BS_aug/BS_ref_aug)
  
  metric_summary[counter,1] <- l
  metric_summary[counter,2] <- cor_jun
  metric_summary[counter,3] <- cor_jul
  metric_summary[counter,4] <- cor_aug
  metric_summary[counter,5] <- tpr_jun
  metric_summary[counter,6] <- tpr_jul
  metric_summary[counter,7] <- tpr_aug
  metric_summary[counter,8] <- tnr_jun
  metric_summary[counter,9] <- tnr_jul
  metric_summary[counter,10] <- tnr_aug
  metric_summary[counter,11] <- acc_jun
  metric_summary[counter,12] <- acc_jul
  metric_summary[counter,13] <- acc_aug
  metric_summary[counter,14] <- v_jun$SEDI
  metric_summary[counter,15] <- v_jul$SEDI
  metric_summary[counter,16] <- v_aug$SEDI
  metric_summary[counter,17] <- BSS_jun
  metric_summary[counter,18] <- BSS_jul
  metric_summary[counter,19] <- BSS_aug
  counter=counter+1
  
}

metric_summary2 <- metric_summary %>% pivot_longer(c(2:4), names_to = "correlation_months", values_to="correlation_values") %>% 
  pivot_longer(c(2:4), names_to = "tpr_months", values_to="tpr_values") %>% 
  pivot_longer(c(2:4), names_to = "tnr_months", values_to="tnr_values") %>% 
  pivot_longer(c(2:4), names_to = "acc_months", values_to="acc_values") %>% 
  pivot_longer(c(2:4), names_to = "sedi_months", values_to="sedi_values") %>% 
  pivot_longer(c(2:4), names_to = "bss_months", values_to="bss_values")
metric_summary2 <- metric_summary2 %>%  separate(correlation_months, c(NA, "correlation_months"))
metric_summary2 <- metric_summary2 %>%  separate(tpr_months, c(NA, "tpr_months"))
metric_summary2 <- metric_summary2 %>%  separate(tnr_months, c(NA, "tnr_months"))
metric_summary2 <- metric_summary2 %>%  separate(acc_months, c(NA, "acc_months"))
metric_summary2 <- metric_summary2 %>%  separate(sedi_months, c(NA, "sedi_months"))
metric_summary2 <- metric_summary2 %>%  separate(bss_months, c(NA, "bss_months"))

saveRDS(metric_summary2,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_8c/metric_summary_global_1980-2010_only3_scenario.rds')

#Correlation plots
metric_summary2$correlation_months <- as.factor(metric_summary2$correlation_months)
metric_summary2$correlation_months <- factor(metric_summary2$correlation_months, levels=rev(levels(metric_summary2$correlation_months)))
cor_plot <- ggplot(data = metric_summary2, aes(x=correlation_months,y=lead,fill=correlation_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(correlation_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.36, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  ggtitle("Correlation")

#TPR plots
metric_summary2$tpr_months <- as.factor(metric_summary2$tpr_months)
metric_summary2$tpr_months <- factor(metric_summary2$tpr_months, levels=rev(levels(metric_summary2$tpr_months)))
tpr_plot <- ggplot(data = metric_summary2, aes(x=tpr_months,y=lead,fill=tpr_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tpr_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TPR") +
  ggtitle("True Positive Rate")

#TNR
metric_summary2$tnr_months <- as.factor(metric_summary2$tnr_months)
metric_summary2$tnr_months <- factor(metric_summary2$tnr_months, levels=rev(levels(metric_summary2$tnr_months)))
tnr_plot <- ggplot(data = metric_summary2, aes(x=tnr_months,y=lead,fill=tnr_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tnr_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TNR") +
  ggtitle("True Negative Rate")

#ACC
metric_summary2$acc_months <- as.factor(metric_summary2$acc_months)
metric_summary2$acc_months <- factor(metric_summary2$acc_months, levels=rev(levels(metric_summary2$acc_months)))
acc_plot <- ggplot(data = metric_summary2, aes(x=acc_months,y=lead,fill=acc_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(acc_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.71, limit = c(0,1), space = "Lab",
                       name="Accuracy") +
  ggtitle("Accuracy")
# (0.18*0.18) + (0.82*0.82) = 0.7048

#SEDI
metric_summary2$sedi_months <- as.factor(metric_summary2$sedi_months)
metric_summary2$sedi_months <- factor(metric_summary2$sedi_months, levels=rev(levels(metric_summary2$sedi_months)))
sedi_plot <- ggplot(data = metric_summary2, aes(x=sedi_months,y=lead,fill=sedi_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(sedi_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="SEDI") +
  ggtitle("SEDI")

#BSS
metric_summary2$bss_months <- as.factor(metric_summary2$bss_months)
metric_summary2$bss_months <- factor(metric_summary2$bss_months, levels=rev(levels(metric_summary2$bss_months)))
bss_plot <- ggplot(data = metric_summary2, aes(x=bss_months,y=lead,fill=bss_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(bss_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="BSS") +
  ggtitle("Brier Skill Score")

tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/fig_skill_plots_global_only3_step8c.tiff', res=200, units="in",width=12,height=8)
grid.arrange(cor_plot, tpr_plot, tnr_plot, acc_plot, sedi_plot, bss_plot, nrow=2)
dev.off()

#-----Step 9b: Estimate Skill: downscaled apples------
#TNR, TPR, ACC, SEDI, and Cor
df_metric <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_7b/total_metric_globalforecasts_1980-2010_apples_ensemble.rds')

metric_summary <- as.data.frame(matrix(NA,nrow=4,ncol=19))
colnames(metric_summary) <- c("lead",
                              "cor_jun","cor_jul","cor_aug",
                              "tpr_jun","tpr_jul","tpr_aug",
                              "tnr_jun","tnr_jul","tnr_aug",
                              "acc_jun","acc_jul","acc_aug",
                              "sedi_jun","sedi_jul","sedi_aug",
                              "bss_jun","bss_jul","bss_aug")

q=unique(df_metric$observed_quantile)
counter=1
for (l in 0.5:5.5){ #start lead time, which then adds 6 months to get lead end time. 
  
  #Correlation coefficients
  cor_jun <- cor(df_metric$june_fcast[df_metric$lead==l], df_metric$june_obs[df_metric$lead==l], use="na.or.complete")
  cor_jul <- cor(df_metric$jul_fcast[df_metric$lead==l], df_metric$jul_obs[ df_metric$lead==l], use="na.or.complete")
  cor_aug <- cor(df_metric$aug_fcast[ df_metric$lead==l], df_metric$aug_obs[ df_metric$lead==l], use="na.or.complete")
  
  #TPR: TP/TP+FN (not x has to be observed)
  tpr_jun <- tpr(o = df_metric$june_obs[ df_metric$lead==l], f = df_metric$june_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tpr_jul <- tpr(o = df_metric$jul_obs[ df_metric$lead==l], f = df_metric$jul_fcast[df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tpr_aug <- tpr(o = df_metric$aug_obs[df_metric$lead==l], f = df_metric$aug_fcast[df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  
  #TNR: TN/TN+FP
  tnr_jun <- tnr(o = df_metric$june_obs[ df_metric$lead==l], f = df_metric$june_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tnr_jul <- tnr(o = df_metric$jul_obs[df_metric$lead==l], f = df_metric$jul_fcast[df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  tnr_aug <- tnr(o = df_metric$aug_obs[ df_metric$lead==l], f = df_metric$aug_fcast [df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  
  #ACC: TN + TN / total
  acc_jun <- acc(o = df_metric$june_obs[ df_metric$lead==l], f = df_metric$june_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  acc_jul <- acc(o = df_metric$jul_obs[df_metric$lead==l], f = df_metric$jul_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  acc_aug <- acc(o = df_metric$aug_obs[ df_metric$lead==l], f = df_metric$aug_fcast[ df_metric$lead==l], qo = q, qf = df_metric$ssta_quantile[df_metric$lead==l])
  
  obs_jun <- df_metric$june_obs[ df_metric$lead==l]
  pred_jun <- df_metric$june_fcast[ df_metric$lead==l]
  obs_jun <- obs_jun[!is.na(obs_jun)]
  pred_jun <- pred_jun[!is.na(pred_jun)]
  v_jun <- verify(obs = ifelse(obs_jun>q,1,0), pred = ifelse(pred_jun>unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0), frcst.type = "binary", obs.type = "binary")
  
  obs_jul <- df_metric$jul_obs[ df_metric$lead==l]
  pred_jul <- df_metric$jul_fcast[df_metric$lead==l]
  obs_jul <- obs_jul[!is.na(obs_jul)]
  pred_jul <- pred_jul[!is.na(pred_jul)]
  v_jul <- verify(obs = ifelse(obs_jul>q,1,0), pred = ifelse(pred_jul>unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0), frcst.type = "binary", obs.type = "binary")
  
  obs_aug <- df_metric$aug_obs[df_metric$lead==l]
  pred_aug <- df_metric$aug_fcast[ df_metric$lead==l]
  obs_aug <- obs_aug[!is.na(obs_aug)]
  pred_aug <- pred_aug[!is.na(pred_aug)]
  v_aug <- verify(obs = ifelse(obs_aug>q,1,0), pred = ifelse(pred_aug>unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0), frcst.type = "binary", obs.type = "binary")
  
  #BS and BSS
  obs_jun <- ifelse(obs_jun <= q,1,0)
  pred_jun <- ifelse(pred_jun <= unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0)
  BS_jun <- mean((pred_jun - obs_jun)^2)
  BS_ref_jun <-  0.82
  BSS_jun <- 1 - (BS_jun/BS_ref_jun)
  
  obs_jul <- ifelse(obs_jul <= q,1,0)
  pred_jul <- ifelse(pred_jul <= unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0)
  BS_jul <- mean((pred_jul - obs_jul)^2)
  BS_ref_jul <-  0.82
  BSS_jul <- 1 - (BS_jul/BS_ref_jul)
  
  obs_aug <- ifelse(obs_aug <= q,1,0)
  pred_aug <- ifelse(pred_aug <= unique(df_metric$ssta_quantile[df_metric$lead==l]),1,0)
  BS_aug <- mean((pred_aug - obs_aug)^2)
  BS_ref_aug <-  0.82
  BSS_aug <- 1 - (BS_aug/BS_ref_aug)
  
  
  metric_summary[counter,1] <- l
  metric_summary[counter,2] <- cor_jun
  metric_summary[counter,3] <- cor_jul
  metric_summary[counter,4] <- cor_aug
  metric_summary[counter,5] <- tpr_jun
  metric_summary[counter,6] <- tpr_jul
  metric_summary[counter,7] <- tpr_aug
  metric_summary[counter,8] <- tnr_jun
  metric_summary[counter,9] <- tnr_jul
  metric_summary[counter,10] <- tnr_aug
  metric_summary[counter,11] <- acc_jun
  metric_summary[counter,12] <- acc_jul
  metric_summary[counter,13] <- acc_aug
  metric_summary[counter,14] <- v_jun$SEDI
  metric_summary[counter,15] <- v_jul$SEDI
  metric_summary[counter,16] <- v_aug$SEDI
  metric_summary[counter,17] <- BSS_jun
  metric_summary[counter,18] <- BSS_jul
  metric_summary[counter,19] <- BSS_aug
  counter=counter+1
  
}

metric_summary2 <- metric_summary %>% pivot_longer(c(2:4), names_to = "correlation_months", values_to="correlation_values") %>% 
  pivot_longer(c(2:4), names_to = "tpr_months", values_to="tpr_values") %>% 
  pivot_longer(c(2:4), names_to = "tnr_months", values_to="tnr_values") %>% 
  pivot_longer(c(2:4), names_to = "acc_months", values_to="acc_values") %>% 
  pivot_longer(c(2:4), names_to = "sedi_months", values_to="sedi_values") %>% 
  pivot_longer(c(2:4), names_to = "bss_months", values_to="bss_values")
metric_summary2 <- metric_summary2 %>%  separate(correlation_months, c(NA, "correlation_months"))
metric_summary2 <- metric_summary2 %>%  separate(tpr_months, c(NA, "tpr_months"))
metric_summary2 <- metric_summary2 %>%  separate(tnr_months, c(NA, "tnr_months"))
metric_summary2 <- metric_summary2 %>%  separate(acc_months, c(NA, "acc_months"))
metric_summary2 <- metric_summary2 %>%  separate(sedi_months, c(NA, "sedi_months"))
metric_summary2 <- metric_summary2 %>%  separate(bss_months, c(NA, "bss_months"))

saveRDS(metric_summary2, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_9b/metric_downscaled.rds')

#Correlation plots
metric_summary2$correlation_months <- as.factor(metric_summary2$correlation_months)
metric_summary2$correlation_months <- factor(metric_summary2$correlation_months, levels=rev(levels(metric_summary2$correlation_months)))
gcor <- ggplot(data = metric_summary2, aes(x=correlation_months,y=lead,fill=correlation_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(correlation_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.36, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  ggtitle("Correlation")

#TPR plots
metric_summary2$tpr_months <- as.factor(metric_summary2$tpr_months)
metric_summary2$tpr_months <- factor(metric_summary2$tpr_months, levels=rev(levels(metric_summary2$tpr_months)))
gtpr <- ggplot(data = metric_summary2, aes(x=tpr_months,y=lead,fill=tpr_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tpr_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TPR") +
  ggtitle("True Positive Rate")

#TNR
metric_summary2$tnr_months <- as.factor(metric_summary2$tnr_months)
metric_summary2$tnr_months <- factor(metric_summary2$tnr_months, levels=rev(levels(metric_summary2$tnr_months)))
gtnr <- ggplot(data = metric_summary2, aes(x=tnr_months,y=lead,fill=tnr_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tnr_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TNR") +
  ggtitle("True Negative Rate")

#ACC
metric_summary2$acc_months <- as.factor(metric_summary2$acc_months)
metric_summary2$acc_months <- factor(metric_summary2$acc_months, levels=rev(levels(metric_summary2$acc_months)))
gacc <- ggplot(data = metric_summary2, aes(x=acc_months,y=lead,fill=acc_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(acc_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.71, limit = c(0,1), space = "Lab",
                       name="Accuracy") +
  ggtitle("Accuracy")
# (0.18*0.18) + (0.82*0.82) = 0.7048

#SEDI
metric_summary2$sedi_months <- as.factor(metric_summary2$sedi_months)
metric_summary2$sedi_months <- factor(metric_summary2$sedi_months, levels=rev(levels(metric_summary2$sedi_months)))
gsedi <- ggplot(data = metric_summary2, aes(x=sedi_months,y=lead,fill=sedi_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(sedi_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="SEDI") +
  ggtitle("SEDI")

#BSS
metric_summary2$bss_months <- as.factor(metric_summary2$bss_months)
metric_summary2$bss_months <- factor(metric_summary2$bss_months, levels=rev(levels(metric_summary2$bss_months)))
gbss <- ggplot(data = metric_summary2, aes(x=bss_months,y=lead,fill=bss_values))+
  geom_tile(color="white")+
  geom_text(aes(label = round(bss_values, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="BSS") +
  ggtitle("Brier Skill Score")

tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/fig_skill_plots_downscaled_apples_step9b.tiff', res=200, units="in",width=12,height=8)
grid.arrange(gcor, gtpr, gtnr, gacc, gsedi,gbss, nrow=2)
dev.off()

#------Step 10a: Compare apples-to-apples in skill----
gl <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_8b/metric_summary_global_1980-2010_applesscenario.rds')
do <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_9b/metric_downscaled.rds')

do$model <- "downscaled"
gl$model <- "global"
all <- rbind(gl, do)

cor_long <- all[,c(1,2,3,14)] %>% distinct() %>% pivot_wider(names_from = model,values_from=correlation_values)
tpr_long <- all[,c(1,4,5,14)]  %>% distinct() %>% pivot_wider(names_from = model,values_from=tpr_values)
tnr_long <- all[,c(1,6,7,14)] %>% distinct() %>% pivot_wider(names_from = model,values_from=tnr_values)
acc_long <- all[,c(1,8,9,14)]  %>% distinct() %>% pivot_wider(names_from = model,values_from=acc_values)
sedi_long <- all[,c(1,10,11,14)]  %>% distinct() %>% pivot_wider(names_from = model,values_from=sedi_values)
bss_long <- all[,c(1,12,13,14)]  %>% distinct() %>% pivot_wider(names_from = model,values_from=bss_values)

cor_long$diff <- cor_long$global-cor_long$downscaled
tpr_long$diff <- tpr_long$global-tpr_long$downscaled
tnr_long$diff <- tnr_long$global-tnr_long$downscaled
acc_long$diff <- acc_long$global-acc_long$downscaled
sedi_long$diff <- sedi_long$global-sedi_long$downscaled
bss_long$diff <- bss_long$global-bss_long$downscaled

cor_long$correlation_months <- as.factor(cor_long$correlation_months)
cor_long$correlation_months <- factor(cor_long$correlation_months, levels=rev(levels(cor_long$correlation_months)))
g1 <- ggplot(cor_long, aes(x=correlation_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.2,0.2), space = "Lab", 
                       name="") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("Differece in Correlation: global - downscaled")

tpr_long$tpr_months <- as.factor(tpr_long$tpr_months)
tpr_long$tpr_months <- factor(tpr_long$tpr_months, levels=rev(levels(tpr_long$tpr_months)))
g2 <- ggplot(tpr_long, aes(x=tpr_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.3,0.3), space = "Lab",
                       name="TPR") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("Differece in TNR: global - downscaled")

tnr_long$tnr_months <- as.factor(tnr_long$tnr_months)
tnr_long$tnr_months <- factor(tnr_long$tnr_months, levels=rev(levels(tnr_long$tnr_months)))
g3 <- ggplot(tnr_long, aes(x=tnr_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.06,0.06), space = "Lab",
                       name="TNR") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("Differece in TPR: global - downscaled")

acc_long$acc_months <- as.factor(acc_long$acc_months)
acc_long$acc_months <- factor(acc_long$acc_months, levels=rev(levels(acc_long$acc_months)))
g4 <- ggplot(acc_long, aes(x=acc_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.05,0.05), space = "Lab",
                       name="Accuracy") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("Differece in ACC: global - downscaled")

sedi_long$sedi_months <- as.factor(sedi_long$sedi_months)
sedi_long$sedi_months <- factor(sedi_long$sedi_months, levels=rev(levels(sedi_long$sedi_months)))
g5 <- ggplot(sedi_long, aes(x=sedi_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.3,0.3), space = "Lab",
                       name="SEDI") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("Differece in SEDI: global - downscaled")

bss_long$bss_months <- as.factor(bss_long$bss_months)
bss_long$bss_months <- factor(bss_long$bss_months, levels=rev(levels(bss_long$bss_months)))
g6 <- ggplot(bss_long, aes(x=bss_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.05,0.05), space = "Lab",
                       name="BSS") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("Differece in BSS: global - downscaled")


tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/fig_skill_plots_comparison_step10.tiff', res=200, units="in",width=15,height=8)
grid.arrange(g1, g3, g2, g4, g5,g6, nrow=2)
dev.off()


#------Step 10b: Compare only3 in skill----
gl <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_8c/metric_summary_global_1980-2010_only3_scenario.rds')
do <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_9b/metric_downscaled.rds')

do$model <- "downscaled"
gl$model <- "global"
all <- rbind(gl, do)

cor_long <- all[,c(1,2,3,14)] %>% distinct() %>% pivot_wider(names_from = model,values_from=correlation_values)
tpr_long <- all[,c(1,4,5,14)]  %>% distinct() %>% pivot_wider(names_from = model,values_from=tpr_values)
tnr_long <- all[,c(1,6,7,14)] %>% distinct() %>% pivot_wider(names_from = model,values_from=tnr_values)
acc_long <- all[,c(1,8,9,14)]  %>% distinct() %>% pivot_wider(names_from = model,values_from=acc_values)
sedi_long <- all[,c(1,10,11,14)]  %>% distinct() %>% pivot_wider(names_from = model,values_from=sedi_values)
bss_long <- all[,c(1,12,13,14)]  %>% distinct() %>% pivot_wider(names_from = model,values_from=bss_values)

cor_long$diff <- cor_long$global-cor_long$downscaled
tpr_long$diff <- tpr_long$global-tpr_long$downscaled
tnr_long$diff <- tnr_long$global-tnr_long$downscaled
acc_long$diff <- acc_long$global-acc_long$downscaled
sedi_long$diff <- sedi_long$global-sedi_long$downscaled
bss_long$diff <- bss_long$global-bss_long$downscaled

cor_long$correlation_months <- as.factor(cor_long$correlation_months)
cor_long$correlation_months <- factor(cor_long$correlation_months, levels=rev(levels(cor_long$correlation_months)))
g1 <- ggplot(cor_long, aes(x=correlation_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.5,0.2), space = "Lab", 
                       name="") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("Differece in Correlation: global - downscaled")

tpr_long$tpr_months <- as.factor(tpr_long$tpr_months)
tpr_long$tpr_months <- factor(tpr_long$tpr_months, levels=rev(levels(tpr_long$tpr_months)))
g2 <- ggplot(tpr_long, aes(x=tpr_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.3,0.3), space = "Lab",
                       name="TPR") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("Differece in TNR: global - downscaled")

tnr_long$tnr_months <- as.factor(tnr_long$tnr_months)
tnr_long$tnr_months <- factor(tnr_long$tnr_months, levels=rev(levels(tnr_long$tnr_months)))
g3 <- ggplot(tnr_long, aes(x=tnr_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.3,0.3), space = "Lab",
                       name="TNR") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("Differece in TPR: global - downscaled")

acc_long$acc_months <- as.factor(acc_long$acc_months)
acc_long$acc_months <- factor(acc_long$acc_months, levels=rev(levels(acc_long$acc_months)))
g4 <- ggplot(acc_long, aes(x=acc_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.15,0.1), space = "Lab",
                       name="Accuracy") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("Differece in ACC: global - downscaled")

sedi_long$sedi_months <- as.factor(sedi_long$sedi_months)
sedi_long$sedi_months <- factor(sedi_long$sedi_months, levels=rev(levels(sedi_long$sedi_months)))
g5 <- ggplot(sedi_long, aes(x=sedi_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.5,0.3), space = "Lab",
                       name="SEDI") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("Differece in SEDI: global - downscaled")

bss_long$bss_months <- as.factor(bss_long$bss_months)
bss_long$bss_months <- factor(bss_long$bss_months, levels=rev(levels(bss_long$bss_months)))
g6 <- ggplot(bss_long, aes(x=bss_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.17,0.05), space = "Lab",
                       name="BSS") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  ggtitle("Differece in BSS: global - downscaled")


tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/fig_skill_plots_comparison_step10b.tiff', res=200, units="in",width=15,height=8)
grid.arrange(g1, g3, g2, g4, g5,g6, nrow=2)
dev.off()


#------Step 11: Case study examples of whether we could forecast-------
# August 2014; June-August 2015; June-August 2016
#Best case scenario from global forecasts

data <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_5a/total_metric_globalforecasts_1980-2020_bestcase_ensemble.rds')

data[data$year==2014,] #leads 0.5, 1.5, 3.5 were correct for August closure
data[data$year==2015,] #ALL leads were correct for June, July, and August closures
data[data$year==2016,] #ALL leads were correct for June, July, and August closures



#----END-------


