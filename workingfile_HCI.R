#Code used to generate forecasts of HCI
#Code used to generate forecasts of HCI
#For full methods and details see: Brodie et al., (in review)
#Data is deposited on dryad: (INSERT LINK)
#Steps below correspond to files stored on dryad

#----librarys-----
library(raster)
library(ncdf4)
library(tidyverse)
library(reshape2)
library(glue)
library(lubridate)
library(maps)
library(mapdata)
library(ggpubr)
library(verification)
library(gridExtra)

#------Step 1a: calculate monthly OBSERVED thresholds------
nc_hist <- nc_open('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/wcra31_sst_monthly_1981_2010.nc')
template=raster("~/Dropbox/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/template.grd")

#Get vars from netcdf
lat <- ncvar_get(nc_hist,'lat')[1,]
lon <- ncvar_get(nc_hist,'lon')[,1]
year_hist <- ncvar_get(nc_hist,'year')
month_hist <- ncvar_get(nc_hist,'month')

#Find dims for HCI area
# Box: -120.3 - -116; 30.8-34.5
lat_dimS <- which(lat==35)
lat_dimN <- which(lat==40)
lon_dimW <- which(lon==-125)#rough, needs to be 75km from shore
lon_dimE <- which(lon==-120)#rough, needs to be 75km from shore

#Get sst
dat_hist <- ncvar_get(nc_hist, 'sst')

#Get monthly mean SST in domain for historical period (1981-2010) 
obs_sst_hist <- as.data.frame(matrix(NA, 360, 3))
colnames(obs_sst_hist) <- c("year","month","sst")
obs_sst_hist$year <- c(year_hist)
obs_sst_hist$month <- c(month_hist)
counter = 1
for (y in 1981:2010){
  for (m in 1:12){
    print(glue("running year {y} and month {m}"))
    
    #Template to force netcdf to raster dims
    template <- raster(xmn=-125,xmx=-120,ymn=35,ymx=40,resolution=0.1)
    
    #Get SST in raster form
    idx_hist <- which(year_hist==y & month_hist==m)
    sst_hist <- dat_hist[lon_dimW:lon_dimE,lat_dimS:lat_dimN, idx_hist]
    lat_dim <- lat[which(lat>=35 & lat<=40)]
    lon_dim <- lon[which(lon>=-125 & lon<=-120)]
    colnames(sst_hist) <- lat_dim
    rownames(sst_hist) <- lon_dim
    sst_hist <- melt(sst_hist)
    sst_r <- rasterFromXYZ(sst_hist)
    sst_r=raster::resample(sst_r,template)
    # plot(sst_r)
    # map('worldHires', add=T, fill=T)
    
    #Get distacne to shore
    d <- readRDS('~/Dropbox/PROJECTS/WRAP Location/r_distcoast.rds')
    d[d[] >= 75] = NA
    d <- crop(d,extent(-125,-120,35,40))
    # plot(d)
    # map('worldHires', add=T, fill=T)
    
    #Now crop sst to distance 
    fin <- mask(sst_r, d)
    # plot(fin)
    # map('worldHires', add=T, fill=T)
    # 
    out <- mean(rasterToPoints(fin)[,3])
    
    obs_sst_hist[counter,1] <- y
    obs_sst_hist[counter,2] <- m
    obs_sst_hist[counter,3] <- out
    
    counter=counter+1
  }
}

head(obs_sst_hist)
thresholds <- obs_sst_hist %>% group_by(month) %>% summarise_at("sst",mean)
# thresholds$sst <- round(thresholds$sst,2)
thresholds$sst <- round(thresholds$sst/0.5) * 0.5
saveRDS(thresholds, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_1a/observed_monthly_thresholds_1980-2010.rds')



#------Step 1b: calculate OBSERVED HCI-------
nc_hist <- nc_open('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/wcra31_sst_monthly_1981_2010.nc')
nc_nrt <- nc_open('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/wcnrt_sst_monthly_201101_202103.nc')
thresholds <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_1a/observed_monthly_thresholds_1980-2010.rds')

#Get vars from netcdf
lat <- ncvar_get(nc_hist,'lat')[1,]
lon <- ncvar_get(nc_hist,'lon')[,1]
year_hist <- ncvar_get(nc_hist,'year')
month_hist <- ncvar_get(nc_hist,'month')
year_nrt <- ncvar_get(nc_nrt,'year')
month_nrt <- ncvar_get(nc_nrt,'month')

#Find dims for HCI area
# Box: -120.3 - -116; 30.8-34.5
lat_dimS <- which(lat==35)
lat_dimN <- which(lat==40)
lon_dimW <- which(lon==-126.5)#rough, needs to be 75km from shore
lon_dimE <- which(lon==-120)#rough, needs to be 75km from shore

#Get sst
dat_hist <- ncvar_get(nc_hist, 'sst')
dat_nrt <- ncvar_get(nc_nrt, 'sst')

#make sure I have data from step 1a
head(obs_sst_hist)
head(thresholds)

#Get monthly mean SST in domain for historical period (1981-2010) 
sst_area <- as.data.frame(matrix(NA, 504, 3))
colnames(sst_area) <- c("year","month","area")
counter = 1
for (y in 1981:2021){
  for (m in 1:12){
    print(glue("running year {y} and month {m}"))
    
    #Template to force netcdf to raster dims
    template <- raster(xmn=-126.5,xmx=-120,ymn=35,ymx=40,resolution=0.1)
    
    if(y <=2010){
    #Get SST in raster form
    idx_hist <- which(year_hist==y & month_hist==m)
    sst_hist <- dat_hist[lon_dimW:lon_dimE,lat_dimS:lat_dimN, idx_hist]
    lat_dim <- lat[which(lat>=35 & lat<=40)]
    lon_dim <- lon[which(lon>=-126.5 & lon<=-120)]
    colnames(sst_hist) <- lat_dim
    rownames(sst_hist) <- lon_dim
    sst_hist <- melt(sst_hist)
    sst_r <- rasterFromXYZ(sst_hist)
    sst_r=raster::resample(sst_r,template)
    # plot(sst_r)
    # map('worldHires', add=T, fill=T)
    
    #Get distacne to shore
    d <- readRDS('~/Dropbox/PROJECTS/WRAP Location/r_distcoast.rds')
    d[d[] >= 150] = NA
    d <- crop(d,extent(-126.5,-120,35,40))
    # plot(d)
    # map('worldHires', add=T, fill=T)
    
    #Now crop sst to distance 
    # fin <- mask(sst_r, d)
    # plot(fin)
    # map('worldHires', add=T, fill=T)
    
    area <- as.data.frame(rasterToPoints(fin))
    out <- nrow(area[area$value <= thresholds$sst[thresholds$month==m],])
    } else {
      #Get SST in raster form
      idx_nrt <- which(year_nrt==y & month_nrt==m)
      sst_nrt <- dat_nrt[lon_dimW:lon_dimE,lat_dimS:lat_dimN, idx_nrt]
      lat_dim <- lat[which(lat>=35 & lat<=40)]
      lon_dim <- lon[which(lon>=-126.5 & lon<=-120)]
      colnames(sst_nrt) <- lat_dim
      rownames(sst_nrt) <- lon_dim
      sst_nrt <- melt(sst_nrt)
      sst_r <- rasterFromXYZ(sst_nrt)
      sst_r=raster::resample(sst_r,template)
      # plot(sst_r)
      # map('worldHires', add=T, fill=T)
      
      #Get distacne to shore
      d <- readRDS('~/Dropbox/PROJECTS/WRAP Location/r_distcoast.rds')
      d[d[] >= 150] = NA
      d <- crop(d,extent(-126.5,-120,35,40))
      # plot(d)
      # map('worldHires', add=T, fill=T)
      
      #Now crop sst to distance 
      fin <- mask(sst_r, d)
      # plot(fin)
      # map('worldHires', add=T, fill=T)
      
      area <- as.data.frame(rasterToPoints(fin))
      out <- nrow(area[area$value <= thresholds$sst[thresholds$month==m],])
    }
    
    
    sst_area[counter,1] <- y
    sst_area[counter,2] <- m
    sst_area[counter,3] <- out
    
    counter=counter+1
  }
}

head(sst_area)

sst_area$area <- sst_area$area/1033
sst_area$date <- as.Date(glue("{sst_area$year}-{sst_area$month}-01"))
mean(sst_area$area, na.rm=T)
sd(sst_area$area, na.rm=T)
plot(sst_area$date,sst_area$area, type="l")
abline(h=0.4, col="blue")
abline(h=0.4+0.28, col="red")
abline(h=0.5-0.28, col="red")

saveRDS(sst_area,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/observed_HCI_1981-2021.rds')

g1 <- ggplot(sst_area, aes(x=date,y=area))+
  geom_line()+
  geom_hline(yintercept = 0.4, linetype="dotted")+
  geom_hline(yintercept = 0.4-0.28, col="blue")+
  geom_hline(yintercept = 0.4+0.28, col="blue")
tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/observed_HCI_timeseries_1981-2021.tiff', units="in", res=300, width=10, height=2)
plot(g1)
dev.off()

#------Step 2a: calculate monthly FORECAST thresholds-----
#Use global ensemble
#CAREFUL: this takes many hours

#First get SST within approximate 150km-from coast box
threshold_function_globalmodels <- function(x,d){
  sst_box <- as.data.frame(matrix(NA, nrow=60000, ncol=6))
  colnames(sst_box) <- c("glob_model","ensemble","init_month","lead_month","forecast_date","sst_mean")
  sst_box$forecast_date <- as.Date(sst_box$forecast_date,"%Y-%m-%d")
  sst_box$init_month <- as.Date(sst_box$init_month,"%Y-%m-%d")
  counter=1
    nc <- x
    nc <- nc_open(nc)
    
    #Get vars from netcdf
    lat <- ncvar_get(nc,'lat')
    lon <- ncvar_get(nc,'lon')
    lead <- ncvar_get(nc,'lead')
    time <- ncvar_get(nc,'time') #months since 1960-01-01. I don't know how to convert this, so next few lines are a workaround
    t1 <- ymd(d)
    time <- time-time[1]
    init_time <- t1 %m+% months(time)
    member <- ncvar_get(nc,'member')
    t_start <- which(init_time>=d)[1]
    t_end <- time[length(time)]
    
    #Find dims
    lat_dimS <- which(lat==35) #has to be 35 not 35.5
    lat_dimN <- which(lat==40)
    lon_dimW <- which(lon==233) #rough, needs to be 75km from shore
    lon_dimE <- which(lon==240) #rough, needs to be 75km from shore
    
    dat <- ncvar_get(nc, 'sst')
    
    for (i in t_start:t_end){ #years/month index
      # for (i in 1:2){ #years/month index
      for (l in 1:length(lead)){ #index
        for (m in 1:length(member)){ #ensemble member index
          
          #Template to force netcdf to raster dims
          template <- raster(xmn=-126.5,xmx=-119.5,ymn=34.5,ymx=40.5,resolution=1)
          
          #Get SST in raster form
          sst_hist <- dat[lon_dimW:lon_dimE,lat_dimS:lat_dimN,m,l,i]
          lat_dim <- lat[which(lat>=35 & lat<=40)]
          lon_dim <- lon[which(lon>=233 & lon<=240)]
          lon_dim <- lon_dim - 360
          colnames(sst_hist) <- lat_dim
          rownames(sst_hist) <- lon_dim
          sst_hist <- melt(sst_hist)
          sst_r <- rasterFromXYZ(sst_hist)
          sst_r=raster::resample(sst_r,template)
          # plot(sst_r)
          # map('worldHires', add=T, fill=T)
          
          #Get distacne to shore. Doesn't quite work!
          # d <- readRDS('~/Dropbox/PROJECTS/WRAP Location/r_distcoast.rds')
          # d[d[] >= 75] = NA
          # d <- crop(d,extent(-125.5,-119.5,34.5,40.5))
          # d=raster::resample(d,template)
          # plot(d,add=T)
          # map('worldHires', add=T, fill=T)
          
          #OR get values at these coordinates:
         coords <- as.data.frame(matrix(NA, nrow = 12,ncol=2))
         colnames(coords) <- c("lon","lat")
         coords$lon <- c(-125,-125,-124,-124,-124,-123,-123,-123,-122,-122,-122,-121)
         coords$lat <- c(40,39,40,39,38,38,37,36,37,36,35,35)
         # points(coords$lon,coords$lat, col="red")
         
          #Now crop sst to distance 
          fin <- raster::extract(sst_r, coords)
    
          out <- mean(fin, na.rm=T)
          out <- round(out/0.5)*0.5
          
          #Get forecast month
          lead.i <- lead[l] - 0.5
          forecast_date <- ymd(init_time[i]) %m+% months(lead.i)
          
          #Write out data
          sst_box[counter,1] <- unlist(strsplit(x,"_"))[7] #global model
          sst_box[counter,2] <- m #ensemble
          sst_box[counter,3] <- init_time[i] #initialisation month index
          sst_box[counter,4] <- lead[l] #lead time 
          sst_box[counter,5] <- forecast_date #get forecast date
          sst_box[counter,6] <- out
          counter=counter+1
        }
      }
    }
  return(sst_box)
}

t1 <- Sys.time()
sst_box_CanCM4 <- threshold_function_globalmodels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_CanCM4i.nc', d= "1981-01-01") #10 ens
saveRDS(sst_box_CanCM4,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/sst_box_CanCM4.rds')
sst_box_COLA <- threshold_function_globalmodels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_COLA-RSMAS-CCSM4.nc', d= "1982-01-01")#10 ens
saveRDS(sst_box_COLA,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/sst_box_COLA.rds')
sst_box_NEMO <- threshold_function_globalmodels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_GEM-NEMO.nc', d= "1981-01-01")#10 ens
saveRDS(sst_box_NEMO,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/sst_box_NEMO.rds')
sst_box_GFDL <- threshold_function_globalmodels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_GFDL-SPEAR.nc', d= "1991-01-01")#15 ens
saveRDS(sst_box_GFDL,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/sst_box_GFDL.rds')
sst_box_NASA <- threshold_function_globalmodels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_NASA-GEOSS2S.nc', d= "1981-02-01") #4 ens
saveRDS(sst_box_NASA,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/sst_box_NASA.rds')
sst_box_NCEP <- threshold_function_globalmodels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_NCEP-CFSv2.nc', d= "1982-01-01") #24 ens
saveRDS(sst_box_NCEP,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/sst_box_NCEP.rds')
t2 <- Sys.time()
t2-t1

sst_box_CanCM4 <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/sst_box_CanCM4.rds')
sst_box_COLA <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/sst_box_COLA.rds')
sst_box_NEMO <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/sst_box_NEMO.rds')
sst_box_GFDL <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/sst_box_GFDL.rds')
sst_box_NASA <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/sst_box_NASA.rds')
sst_box_NCEP <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/sst_box_NCEP.rds')

dat_fcasts <- rbind(sst_box_CanCM4, sst_box_COLA, sst_box_NEMO,
                    sst_box_GFDL, sst_box_NASA, sst_box_NCEP) 
dat_fcasts <- na.omit(dat_fcasts)

#Create threshold based off monthly clim from 1981-2010
dat_fcasts$target_month <- lubridate::month(dat_fcasts$forecast_date)
dat_fcasts$init_month <- lubridate::month(dat_fcasts$init_month)
thresh_forecasts <- dat_fcasts[dat_fcasts$forecast_date<="2010-12-31" & dat_fcasts$forecast_date>="1981-01-01",] %>% 
  group_by(glob_model, ensemble, lead_month, init_month, target_month) %>% summarise_at(vars("sst_mean"),funs(mean(., na.rm=T))) 
colnames(thresh_forecasts)[6] <- "sst_clim"
# thresh_forecasts$sst_clim <- ifelse(thresh_forecasts$sst_clim==0,NA,thresh_forecasts$sst_clim)
ggplot(thresh_forecasts, aes(x=target_month, y=sst_clim, group=ensemble))+
  geom_point(aes(col=glob_model),shape=1)+
  facet_wrap(~lead_month)+
  theme(legend.position = "bottom")

saveRDS(thresh_forecasts, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/forecast_monthly_thresholds_1980-2010.rds')


#------Step 2b: calculate FORECAST HCI-------
#Number of grid cells less than the monthly threshold

thresh_forecasts <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2a/forecast_monthly_thresholds_1980-2010.rds')

HCI_function_globalmodels <- function(x,d){
  fcast_area <- as.data.frame(matrix(NA, nrow=60000, ncol=6))
  colnames(fcast_area) <- c("glob_model","ensemble","init_month","lead_month","forecast_date","sst_area")
  fcast_area$forecast_date <- as.Date(fcast_area$forecast_date,"%Y-%m-%d")
  fcast_area$init_month <- as.Date(fcast_area$init_month,"%Y-%m-%d")
  counter=1
  nc <- x
  nc <- nc_open(nc)
  
  #Get vars from netcdf
  lat <- ncvar_get(nc,'lat')
  lon <- ncvar_get(nc,'lon')
  lead <- ncvar_get(nc,'lead')
  time <- ncvar_get(nc,'time') #months since 1960-01-01. I don't know how to convert this, so next few lines are a workaround
  t1 <- ymd(d)
  time <- time-time[1]
  init_time <- t1 %m+% months(time)
  member <- ncvar_get(nc,'member')
  t_start <- which(init_time>=d)[1]
  t_end <- time[length(time)]
  
  #Find dims
  lat_dimS <- which(lat==35) #has to be 35 not 35.5
  lat_dimN <- which(lat==40)
  lon_dimW <- which(lon==233) #rough, needs to be 150km from shore
  lon_dimE <- which(lon==240) #rough, needs to be 75km from shore
  
  dat <- ncvar_get(nc, 'sst')
  
  for (i in t_start:t_end){ #years/month index
    # for (i in 1:2){ #years/month index
    for (l in 1:length(lead)){ #index
      for (m in 1:length(member)){ #ensemble member index
        
        #Template to force netcdf to raster dims
        template <- raster(xmn=-126.5,xmx=-119.5,ymn=34.5,ymx=40.5,resolution=1)
        
        #Get SST in raster form
        sst_hist <- dat[lon_dimW:lon_dimE,lat_dimS:lat_dimN,m,l,i]
        lat_dim <- lat[which(lat>=35 & lat<=40)]
        lon_dim <- lon[which(lon>=233 & lon<=240)]
        lon_dim <- lon_dim - 360
        colnames(sst_hist) <- lat_dim
        rownames(sst_hist) <- lon_dim
        sst_hist <- melt(sst_hist)
        sst_r <- rasterFromXYZ(sst_hist)
        sst_r=raster::resample(sst_r,template)
        # plot(sst_r)
        # map('worldHires', add=T, fill=T)
        
        #Get distacne to shore
        # d <- readRDS('~/Dropbox/PROJECTS/WRAP Location/r_distcoast.rds')
        # d[d[] >= 150] = NA
        # d <- crop(d,extent(-126.5,-119.5,34.5,40.5))
        # d=raster::resample(d,template)
        # plot(d, add=T)
        # map('worldHires', add=T, fill=T)
        # 
        
        #OR get values at these coordinates:
        coords <- as.data.frame(matrix(NA, nrow = 18,ncol=2))
        colnames(coords) <- c("lon","lat")
        coords$lon <- c(-126,-126,-125,-125,-125,-124,-124,-124,-124,-124,-123,-123,-123,-123,-122,-122,-122,-121)
        coords$lat <- c(40,39,40,39,38,37,40,39,38,36,38,37,36,35,37,36,35,35)
        # points(coords$lon,coords$lat, col="red")
        
        #Now crop sst to distance 
        # fin <- mask(sst_r, d)
        fin <- raster::extract(sst_r, coords)
        # plot(fin)
        # map('worldHires', add=T, fill=T)
        
        #Get forecast month
        lead.i <- lead[l] - 0.5
        forecast_date <- ymd(init_time[i]) %m+% months(lead.i)
        forecast_month <- lubridate::month(forecast_date)
        init_month <- lubridate::month(init_time[i])
        
        
        # area <- as.data.frame(rasterToPoints(fin))
        g <- unlist(strsplit(x,"_"))[7] #global model
        t <- thresh_forecasts$sst_clim[thresh_forecasts$ensemble==m &
                                         thresh_forecasts$lead_month == lead[l] & 
                                         thresh_forecasts$glob_model==g & 
                                         thresh_forecasts$init_month==init_month & 
                                         thresh_forecasts$target_month==forecast_month]
        fin <- na.omit(fin)
        out <- length(fin[fin <= t])
        
        #Write out data
        fcast_area[counter,1] <- g #global model
        fcast_area[counter,2] <- m #ensemble
        fcast_area[counter,3] <- init_time[i] #initialisation month index
        fcast_area[counter,4] <- lead[l] #lead time 
        fcast_area[counter,5] <- forecast_date #get forecast date
        fcast_area[counter,6] <- out
        counter=counter+1
      }
    }
  }
  return(fcast_area)
}

t1 <- Sys.time()
fcast_area_CanCM4 <- HCI_function_globalmodels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_CanCM4i.nc', d= "1981-01-01") #10 ens
saveRDS(fcast_area_CanCM4,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/fcast_area_CanCM4.rds')
fcast_area_COLA <- HCI_function_globalmodels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_COLA-RSMAS-CCSM4.nc', d= "1982-01-01")#10 ens
saveRDS(fcast_area_COLA,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/fcast_area_COLA.rds')
fcast_area_NEMO <- HCI_function_globalmodels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_GEM-NEMO.nc', d= "1981-01-01")#10 ens
saveRDS(fcast_area_NEMO,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/fcast_area_NEMO.rds')
fcast_area_GFDL <- HCI_function_globalmodels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_GFDL-SPEAR.nc', d= "1991-01-01")#15 ens
saveRDS(fcast_area_GFDL,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/fcast_area_GFDL.rds')
fcast_area_NASA <- HCI_function_globalmodels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_NASA-GEOSS2S.nc', d= "1981-02-01") #4 ens
saveRDS(fcast_area_NASA,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/fcast_area_NASA.rds')
fcast_area_NCEP <- HCI_function_globalmodels(x='~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_NCEP-CFSv2.nc', d= "1982-01-01") #24 ens
saveRDS(fcast_area_NCEP,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/fcast_area_NCEP.rds')
t2 <- Sys.time()
t2-t1

summary(fcast_area_CanCM4)
summary(fcast_area_COLA)
summary(fcast_area_NEMO)
summary(fcast_area_GFDL)
summary(fcast_area_NASA)
summary(fcast_area_NCEP)
#GFDL has 8 cells; NCEP has 14; CANCM4 has 13 cells; COLA has 11 cells; NEMO has 13 cells; NASA has 14 cells

fcast_area_CanCM4$area <- fcast_area_CanCM4$sst_area/max(fcast_area_CanCM4$sst_area, na.rm=T)
fcast_area_COLA$area <- fcast_area_COLA$sst_area/max(fcast_area_COLA$sst_area, na.rm=T)
fcast_area_NEMO$area <- fcast_area_NEMO$sst_area/max(fcast_area_NEMO$sst_area, na.rm=T)
fcast_area_GFDL$area <- fcast_area_GFDL$sst_area/max(fcast_area_GFDL$sst_area, na.rm=T)
fcast_area_NASA$area <- fcast_area_NASA$sst_area/max(fcast_area_NASA$sst_area, na.rm=T)
fcast_area_NCEP$area <- fcast_area_NCEP$sst_area/max(fcast_area_NCEP$sst_area, na.rm=T)

fcast_area <- rbind(fcast_area_CanCM4, fcast_area_COLA, fcast_area_GFDL, fcast_area_NASA, fcast_area_NCEP, fcast_area_NEMO)
fcast_area <- na.omit(fcast_area)
saveRDS(fcast_area,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/forecast_allmembers_HCI_1981-2021.rds')

fcast_area_ensemble <- fcast_area %>% group_by(init_month, lead_month, forecast_date) %>% 
                          summarise_at("area", mean)
saveRDS(fcast_area_ensemble,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/forecast_ensmean_HCI_1981-2021.rds')


mean(fcast_area_ensemble$area, na.rm=T)
sd(fcast_area_ensemble$area, na.rm=T)
g1 <- ggplot(data = fcast_area_ensemble, aes(x=forecast_date,y=area))+
  geom_line()+
  geom_hline(yintercept=0.43, linetype="dotted")+
  geom_hline(yintercept=0.43+0.12, col="blue")+
  geom_hline(yintercept=0.43-0.12, col="blue")+
  facet_wrap(~lead_month)
 

tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/forecast_HCI_timeseries_1981-2021.tiff', units="in", res=300, width=12, height=9)
plot(g1)
dev.off()

#------Step 2c: save HCI thresholds comparative scenario -----
fcast_area_ensemble <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/forecast_ensmean_HCI_1981-2021.rds')

me <- fcast_area_ensemble %>% group_by(lead_month) %>% summarise_at("area", mean)
sd <- fcast_area_ensemble %>% group_by(lead_month) %>% summarise_at("area", sd)
thresholds <-  left_join(me, sd, by = "lead_month")
colnames(thresholds) <- c("lead_month", "mean", "sd")
thresholds$high_hci <- thresholds$mean - thresholds$sd
saveRDS(thresholds, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2c/global_forecast_hci_thresholds.rds')

#------Step 2d: save HCI thresholds comparative scenario----
fcast_area <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/forecast_allmembers_HCI_1981-2021.rds')
fcast_area <- fcast_area[fcast_area$glob_model=="CanCM4i.nc",]
fcast_area <- fcast_area[fcast_area$ensemble==2 | fcast_area$ensemble==8 |fcast_area$ensemble==10,]
fcast_area_ensemble_only3 <- fcast_area %>% group_by(init_month, lead_month, forecast_date) %>% 
  summarise_at("area", mean)
saveRDS(fcast_area_ensemble_only3,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/forecast_ensmean_only3_HCI_1981-2010.rds')

me <- fcast_area_ensemble_only3 %>% group_by(lead_month) %>% summarise_at("area", mean)
sd <- fcast_area_ensemble_only3 %>% group_by(lead_month) %>% summarise_at("area", sd)
thresholds <-  left_join(me, sd, by = "lead_month")
colnames(thresholds) <- c("lead_month", "mean", "sd")
thresholds$high_hci <- thresholds$mean - thresholds$sd
saveRDS(thresholds, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2d/global_forecast_hci_thresholds.rds')


#------SKILL FUNCTIONS: tnr, tpr, accuracy----
#Skill Functions

#tp = true positive
#tn = true negative
#fp = false positive (forecast but not observed, i.e. predict high compression but not observed)
#fn = false negative (observed but not forecast, i.e. observed high compression but didn't predict it)

tpr <- function(o,f,qo,qf){ #observed, forecast, observed quantile, and forecast quantile
  tp <- length(which(o<=qo & f<=qf)) #Caution had to flip the logic of these for HCI
  tn <- length(which(o>qo & f>qf))
  fp <- length(which(o>qo & f<qf))
  fn <- length(which(o<qo & f>qf))
  tpr <- tp/(tp+fn)
  return(tpr)
}

tnr <- function(o,f,qo,qf){
  tp <- length(which(o<=qo & f<=qf))#Caution had to flip the logic of these for HCI
  tn <- length(which(o>qo & f>qf))
  fp <- length(which(o>qo & f<qf))
  fn <- length(which(o<qo & f>qf))
  tnr <- tn/(tn+fp)
  return(tnr)
}

acc <- function(o,f,qo,qf){
  tp <- length(which(o<=qo & f<=qf))#Caution had to flip the logic of these for HCI
  tn <- length(which(o>qo & f>qf))
  fp <- length(which(o>qo & f<qf))
  fn <- length(which(o<qo & f>qf))
  acc <- (tp + tn) / sum(tp,tn,fp,fn)
  return(acc)
}





#------Step 3a: evaluate forecast skill: global models all years------
#Correlation; TNR; TPR: Accuracy; SEDI; Brier

#Merge in full data
fcast_area_ensemble <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/forecast_ensmean_HCI_1981-2021.rds')
obs_area <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_1b/observed_HCI_1981-2021.rds')
fcast_area_ensemble$date <- fcast_area_ensemble$forecast_date
colnames(fcast_area_ensemble)[4] <- "fcast_area"
colnames(obs_area)[3] <- "obs_area"
head(fcast_area_ensemble)
head(obs_area)
data <- left_join(fcast_area_ensemble, obs_area[,3:4], "date")
data$month <- lubridate::month(data$date)
head(data)

#Plot correlations
ggplot(data, aes(x=fcast_area, y=obs_area))+
  geom_point()+  
  stat_cor(method = "pearson", label.x = 0, label.y = 1, size=3)+
  facet_wrap(~lead_month)


#Compute all skill metrics
metric_summary <- as.data.frame(matrix(NA,nrow=144,ncol=9))
colnames(metric_summary) <- c("month","lead","cor","cor_sig",
                              "tpr", "tnr","acc","sedi","bss")
counter=1
for (l in 0.5:11.5){ #start lead time, which then adds 6 months to get lead end time. 
  for (m in 1:12){
  
  data_temp <- data[data$lead_month==l & data$month==m,]
  
  #Correlation coefficients
  cor_i <- cor(data_temp$fcast_area, data_temp$obs_area, use="na.or.complete")
  
  #Calculate N effective degrees of freedom to determine if correlation is significant 
  N <- length(unique(lubridate::year(data_temp$forecast_date)))
  tau <- 0:(N-1)
  dm <- list()
  for (t in 1:length(tau)){
    if(t<N-1){
      rx <- cor(data_temp$fcast_area[1:(N-tau[t])],data_temp$fcast_area[(tau[t]+1):N])
      ry <- cor(data_temp$obs_area[1:(N-tau[t])],data_temp$obs_area[(tau[t]+1):N], use="na.or.complete")
    } else {
      rx <- 1
      ry <- 1
    }
    dm[[t]] <- (1-tau[t]/N) * rx * ry
  }
  dm <- unlist(dm)
  NEFF <- N/sum(dm, na.rm=T)  
  # Z <- 0.5 * log((1+cor_i)/(1-cor_i))
  # FisherZ(cor_i) #does the same as above
  # FisherZInv(Z) #just a check
  ci <- CorCI(cor_i,n=NEFF, conf.level = 0.95)
  ci_l <- ifelse(ci[2]>0,1,0 )
  
  #Quantile threshold
  qo = mean(data$obs_area[data$lead_month==l], na.rm=T) #Mean of full time-series
  qf = mean(data$fcast_area[data$lead_month==l], na.rm=T)
  
  #TPR: TP/TP+FN
  tpr_i <- tpr(o = data$obs_area[data$lead_month==l & data$month==m], f = data$fcast_area[data$lead_month==l & data$month==m], qo = qo, qf = qf)
  
  #TNR: TN/TN+FP
  tnr_i <- tnr(o = data$obs_area[data$lead_month==l & data$month==m], f = data$fcast_area[data$lead_month==l & data$month==m], qo = qo, qf = qf)
  
  #ACC: TN + TN / total
  acc_i <- acc(o = data$obs_area[data$lead_month==l & data$month==m], f = data$fcast_area[data$lead_month==l & data$month==m], qo = qo, qf = qf)

  #SEDI
  obs <- data$obs_area[data$lead_month==l & data$month==m]
  pred <-  data$fcast_area[data$lead_month==l & data$month==m]
  if (any(is.na(obs))){
  na_obs <- which(is.na(obs))
  obs <- na.omit(obs)
  pred <- pred[-na_obs]
  }
  v <- verify(obs = ifelse(obs<=qo,1,0), pred = ifelse(pred<=qf,1,0), frcst.type = "binary", obs.type = "binary")

  #BS and BSS
  obs_bi <- ifelse(obs <= qo,1,0)
  pred_bi <- ifelse(pred <= qf,1,0)
  BS <- mean((pred_bi - obs_bi)^2)
  BS_ref <-  0.5
  BSS <- 1 - (BS/BS_ref)
  
  metric_summary[counter,1] <- m
  metric_summary[counter,2] <- l
  metric_summary[counter,3] <- cor_i
  metric_summary[counter,4] <- ci_l
  metric_summary[counter,5] <- tpr_i
  metric_summary[counter,6] <- tnr_i
  metric_summary[counter,7] <- acc_i
  metric_summary[counter,8] <- v$SEDI
  metric_summary[counter,9] <- BSS
  counter=counter+1
  }
}

saveRDS(metric_summary,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_3a/global_forecasts_metric_summary.rds')

#Now plot
#Correlation plots
cor <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=cor))+
  geom_tile(color="white")+
  geom_text(aes(label = cor_sig)) +
  # geom_text(aes(label = round(cor, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson\nCorrelation") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Correlation")

#TPR plots
tpr <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=tpr))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tpr, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TPR") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("TPR")

#TNR
tnr <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=tnr))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tnr, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TNR") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("TNR")

#ACC
acc <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=acc))+
  geom_tile(color="white")+
  geom_text(aes(label = round(acc, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="Accuracy") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Accuracy")
# (0.21*0.21) + (0.79*0.79) = 0.67
# (0.5*0.5) + (0.5*0.5) = 0.5 #based on median (50% quantile) being roughly 0.5

#SEDI
sedi <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=sedi))+
  geom_tile(color="white")+
  geom_text(aes(label = round(sedi, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="SEDI") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("SEDI")

#BSS
bss <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=bss))+
  geom_tile(color="white")+
  geom_text(aes(label = round(bss, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(0,1), space = "Lab",
                       name="BSS") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("BSS")

tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_3a/fig_skill_plots.tiff', res=200, units="in",width=15,height=8)
grid.arrange(cor, tpr, tnr, acc, sedi,bss, nrow=2)
dev.off()

#------Step 3b: evaluate forecast skill: global models 1980-2010 comparative scenario------
#Correlation; TNR; TPR: Accuracy; SEDI; Brier

#Merge in full data
fcast_area_ensemble <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/forecast_ensmean_HCI_1981-2021.rds')
obs_area <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_1b/observed_HCI_1981-2021.rds')
fcast_area_ensemble$date <- fcast_area_ensemble$forecast_date
colnames(fcast_area_ensemble)[4] <- "fcast_area"
colnames(obs_area)[3] <- "obs_area"
head(fcast_area_ensemble)
head(obs_area)
data <- left_join(fcast_area_ensemble, obs_area[,3:4], "date")
data$month <- lubridate::month(data$date)
head(data)
data <- data[data$forecast_date<="2011-12-01",] #Reduce to match downscaled years

#Compute all skill metrics
metric_summary <- as.data.frame(matrix(NA,nrow=144,ncol=9))
colnames(metric_summary) <- c("month","lead","cor","cor_sig",
                              "tpr", "tnr","acc","sedi","bss")
counter=1
for (l in 0.5:11.5){ #start lead time, which then adds 6 months to get lead end time. 
  for (m in 1:12){
    
    data_temp <- data[data$lead_month==l & data$month==m,]
    
    #Correlation coefficients
    cor_i <- cor(data_temp$fcast_area, data_temp$obs_area, use="na.or.complete")
    
    #Calculate N effective degrees of freedom to determine if correlation is significant 
    N <- length(unique(lubridate::year(data_temp$forecast_date)))
    tau <- 0:(N-1)
    dm <- list()
    for (t in 1:length(tau)){
      if(t<N-1){
        rx <- cor(data_temp$fcast_area[1:(N-tau[t])],data_temp$fcast_area[(tau[t]+1):N])
        ry <- cor(data_temp$obs_area[1:(N-tau[t])],data_temp$obs_area[(tau[t]+1):N], use="na.or.complete")
      } else {
        rx <- 1
        ry <- 1
      }
      dm[[t]] <- (1-tau[t]/N) * rx * ry
    }
    dm <- unlist(dm)
    NEFF <- N/sum(dm, na.rm=T)  
    # Z <- 0.5 * log((1+cor_i)/(1-cor_i))
    # FisherZ(cor_i) #does the same as above
    # FisherZInv(Z) #just a check
    ci <- CorCI(cor_i,n=NEFF, conf.level = 0.95)
    ci_l <- ifelse(ci[2]>0,1,0 )
    
    #Quantile threshold
    qo = mean(data$obs_area[data$lead_month==l], na.rm=T) #Mean of full time-series
    qf = mean(data$fcast_area[data$lead_month==l], na.rm=T)
    
    #TPR: TP/TP+FN
    tpr_i <- tpr(o = data$obs_area[data$lead_month==l & data$month==m], f = data$fcast_area[data$lead_month==l & data$month==m], qo = qo, qf = qf)
    
    #TNR: TN/TN+FP
    tnr_i <- tnr(o = data$obs_area[data$lead_month==l & data$month==m], f = data$fcast_area[data$lead_month==l & data$month==m], qo = qo, qf = qf)
    
    #ACC: TN + TN / total
    acc_i <- acc(o = data$obs_area[data$lead_month==l & data$month==m], f = data$fcast_area[data$lead_month==l & data$month==m], qo = qo, qf = qf)
    
    #SEDI
    obs <- data$obs_area[data$lead_month==l & data$month==m]
    pred <-  data$fcast_area[data$lead_month==l & data$month==m]
    if (any(is.na(obs))){
      na_obs <- which(is.na(obs))
      obs <- na.omit(obs)
      pred <- pred[-na_obs]
    }
    v <- verify(obs = ifelse(obs<=qo,1,0), pred = ifelse(pred<=qf,1,0), frcst.type = "binary", obs.type = "binary")
    
    #BS and BSS
    obs_bi <- ifelse(obs <= qo,1,0)
    pred_bi <- ifelse(pred <= qf,1,0)
    BS <- mean((pred_bi - obs_bi)^2)
    BS_ref <-  0.5
    BSS <- 1 - (BS/BS_ref)
    
    metric_summary[counter,1] <- m
    metric_summary[counter,2] <- l
    metric_summary[counter,3] <- cor_i
    metric_summary[counter,4] <- ci_l
    metric_summary[counter,5] <- tpr_i
    metric_summary[counter,6] <- tnr_i
    metric_summary[counter,7] <- acc_i
    metric_summary[counter,8] <- v$SEDI
    metric_summary[counter,9] <- BSS
    counter=counter+1
  }
}

saveRDS(metric_summary,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_3b/global_forecasts_metric_summary.rds')

#Now plot
#Correlation plots
cor <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=cor))+
  geom_tile(color="white")+
  geom_text(aes(label = round(cor, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.38, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Correlation")

#TPR plots
tpr <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=tpr))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tpr, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TPR") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("TPR")

#TNR
tnr <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=tnr))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tnr, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TNR") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("TNR")

#ACC
acc <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=acc))+
  geom_tile(color="white")+
  geom_text(aes(label = round(acc, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="Accuracy") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Accuracy")
# (0.21*0.21) + (0.79*0.79) = 0.67
# (0.5*0.5) + (0.5*0.5) = 0.5 #based on median (50% quantile) being roughly 0.5

#SEDI
sedi <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=sedi))+
  geom_tile(color="white")+
  geom_text(aes(label = round(sedi, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="SEDI") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("SEDI")

#BSS
bss <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=bss))+
  geom_tile(color="white")+
  geom_text(aes(label = round(bss, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(0,1), space = "Lab",
                       name="BSS") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("BSS")

tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_3b/fig_skill_plots.tiff', res=200, units="in",width=15,height=8)
grid.arrange(cor, tpr, tnr, acc, sedi,bss, nrow=2)
dev.off()

#------Step 3c: evaluate forecast skill: global models 1980-2010 reduced scenario------
#Correlation; TNR; TPR: Accuracy; SEDI; Brier

#Merge in full data
fcast_area_ensemble <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/forecast_ensmean_only3_HCI_1981-2010.rds')
obs_area <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_1b/observed_HCI_1981-2021.rds')
fcast_area_ensemble$date <- fcast_area_ensemble$forecast_date
colnames(fcast_area_ensemble)[4] <- "fcast_area"
colnames(obs_area)[3] <- "obs_area"
head(fcast_area_ensemble)
head(obs_area)
data <- left_join(fcast_area_ensemble, obs_area[,3:4], "date")
data$month <- lubridate::month(data$date)
head(data)
data <- data[data$forecast_date<="2011-12-01",] #Reduce to match downscaled years

#Compute all skill metrics
metric_summary <- as.data.frame(matrix(NA,nrow=144,ncol=9))
colnames(metric_summary) <- c("month","lead","cor","cor_sig",
                              "tpr", "tnr","acc","sedi","bss")
counter=1
for (l in 0.5:11.5){ #start lead time, which then adds 6 months to get lead end time. 
  for (m in 1:12){
    
    data_temp <- data[data$lead_month==l & data$month==m,]
    
    #Correlation coefficients
    cor_i <- cor(data_temp$fcast_area, data_temp$obs_area, use="na.or.complete")
    
    #Calculate N effective degrees of freedom to determine if correlation is significant 
    N <- length(unique(lubridate::year(data_temp$forecast_date)))
    tau <- 0:(N-1)
    dm <- list()
    for (t in 1:length(tau)){
      if(t<N-1){
        rx <- cor(data_temp$fcast_area[1:(N-tau[t])],data_temp$fcast_area[(tau[t]+1):N])
        ry <- cor(data_temp$obs_area[1:(N-tau[t])],data_temp$obs_area[(tau[t]+1):N], use="na.or.complete")
      } else {
        rx <- 1
        ry <- 1
      }
      dm[[t]] <- (1-tau[t]/N) * rx * ry
    }
    dm <- unlist(dm)
    NEFF <- N/sum(dm, na.rm=T)  
    # Z <- 0.5 * log((1+cor_i)/(1-cor_i))
    # FisherZ(cor_i) #does the same as above
    # FisherZInv(Z) #just a check
    ci <- CorCI(cor_i,n=NEFF, conf.level = 0.95)
    ci_l <- ifelse(ci[2]>0,1,0 )
    
    #Quantile threshold
    qo = mean(data$obs_area[data$lead_month==l], na.rm=T) #Mean of full time-series
    qf = mean(data$fcast_area[data$lead_month==l], na.rm=T)
    
    #TPR: TP/TP+FN
    tpr_i <- tpr(o = data$obs_area[data$lead_month==l & data$month==m], f = data$fcast_area[data$lead_month==l & data$month==m], qo = qo, qf = qf)
    
    #TNR: TN/TN+FP
    tnr_i <- tnr(o = data$obs_area[data$lead_month==l & data$month==m], f = data$fcast_area[data$lead_month==l & data$month==m], qo = qo, qf = qf)
    
    #ACC: TN + TN / total
    acc_i <- acc(o = data$obs_area[data$lead_month==l & data$month==m], f = data$fcast_area[data$lead_month==l & data$month==m], qo = qo, qf = qf)
    
    #SEDI
    obs <- data$obs_area[data$lead_month==l & data$month==m]
    pred <-  data$fcast_area[data$lead_month==l & data$month==m]
    if (any(is.na(obs))){
      na_obs <- which(is.na(obs))
      obs <- na.omit(obs)
      pred <- pred[-na_obs]
    }
    v <- verify(obs = ifelse(obs<=qo,1,0), pred = ifelse(pred<=qf,1,0), frcst.type = "binary", obs.type = "binary")
    
    #BS and BSS
    obs_bi <- ifelse(obs <= qo,1,0)
    pred_bi <- ifelse(pred <= qf,1,0)
    BS <- mean((pred_bi - obs_bi)^2)
    BS_ref <-  0.5
    BSS <- 1 - (BS/BS_ref)
    
    metric_summary[counter,1] <- m
    metric_summary[counter,2] <- l
    metric_summary[counter,3] <- cor_i
    metric_summary[counter,4] <- ci_l
    metric_summary[counter,5] <- tpr_i
    metric_summary[counter,6] <- tnr_i
    metric_summary[counter,7] <- acc_i
    metric_summary[counter,8] <- v$SEDI
    metric_summary[counter,9] <- BSS
    counter=counter+1
  }
}

saveRDS(metric_summary,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_3c/global_forecasts_metric_summary_only3.rds')

#Now plot
#Correlation plots
corp <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=cor))+
  geom_tile(color="white")+
  geom_text(aes(label = round(cor, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.38, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Correlation")

#TPR plots
tprp <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=tpr))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tpr, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TPR") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("TPR")

#TNR
tnrp <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=tnr))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tnr, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TNR") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("TNR")

#ACC
accp <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=acc))+
  geom_tile(color="white")+
  geom_text(aes(label = round(acc, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="Accuracy") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Accuracy")
# (0.21*0.21) + (0.79*0.79) = 0.67
# (0.5*0.5) + (0.5*0.5) = 0.5 #based on median (50% quantile) being roughly 0.5

#SEDI
sedip <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=sedi))+
  geom_tile(color="white")+
  geom_text(aes(label = round(sedi, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="SEDI") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("SEDI")

#BSS
bssp <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=bss))+
  geom_tile(color="white")+
  geom_text(aes(label = round(bss, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(0,1), space = "Lab",
                       name="BSS") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("BSS")

tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_3c/fig_skill_plots.tiff', res=200, units="in",width=15,height=8)
grid.arrange(corp, tprp, tnrp, accp, sedip,bssp, nrow=2)
dev.off()

#-----Step 4a: calculate monthly DOWNSCALED FORECAST thresholds: -----
#Use downscaled forecasts
#Netcdf of forecasts stores on dropbox
init1_files <- list.files('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/ForecastFields/', full.names = T, pattern="init1")
init7_files <- list.files('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/ForecastFields/', full.names = T, pattern="init7")
sst_files_init1 <- init1_files[grep('sst',init1_files)]
sst_files_init7 <- init7_files[grep('sst',init7_files)]

#First get SST within approximate 150km-from coast box
threshold_function_downscaledmodels <- function(x){
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
    
    #Find dims for HCI area
    # Box: -120.3 - -116; 30.8-34.5
    lat_dimS <- which(lat==35)
    lat_dimN <- which(lat==40)
    lon_dimW <- which(lon==-125)#rough, needs to be 75km from shore
    lon_dimE <- which(lon==-120)#rough, needs to be 75km from shore
    
    #Get monthly climatology from 2003-2010
    dat <- ncvar_get(nc, 'sst')
    
    for (l in 1:12){ #lead times
      for (y in 1:29){ #years (because each forecast initialization on January OR July)
        
        #Template to force netcdf to raster dims
        template <- raster(xmn=-125,xmx=-120,ymn=35,ymx=40,resolution=0.1)
        
        #Get SST in raster form
        sst_hist <- dat[lon_dimW:lon_dimE,lat_dimS:lat_dimN,l,y]
        lat_dim <- lat[which(lat>=35 & lat<=40)]
        lon_dim <- lon[which(lon>=-125 & lon<=-120)]
        colnames(sst_hist) <- lat_dim
        rownames(sst_hist) <- lon_dim
        sst_hist <- melt(sst_hist)
        sst_r <- rasterFromXYZ(sst_hist)
        sst_r=raster::resample(sst_r,template)
        # plot(sst_r)
        # map('worldHires', add=T, fill=T)
        
        #Get distacne to shore
        d <- readRDS('~/Dropbox/PROJECTS/WRAP Location/r_distcoast.rds')
        d[d[] >= 75] = NA
        d <- crop(d,extent(-125,-120,35,40))
        # plot(d)
        # map('worldHires', add=T, fill=T)
        
        #Now crop sst to distance 
        fin <- mask(sst_r, d)
        # plot(fin)
        # map('worldHires', add=T, fill=T)
        # 
        out <- mean(rasterToPoints(fin)[,3])
        
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
        sst_box[counter,5] <- out
        counter=counter+1
      }
    }
  }
  return(sst_box)
}

t1 <- Sys.time()
sst_box_init1 <- threshold_function_downscaledmodels(x=sst_files_init1)
saveRDS(sst_box_init1,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_4a/sst_box_init1.rds')
sst_box_init7 <- threshold_function_downscaledmodels(x=sst_files_init7)
saveRDS(sst_box_init7,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_4a/sst_box_init7.rds')
t2 <- Sys.time()
t2-t1

sst_box_init1 <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_4a/sst_box_init1.rds')
sst_box_init7 <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_4a/sst_box_init7.rds')

dat_fcasts <- rbind(sst_box_init1, sst_box_init7) 
dat_fcasts <- na.omit(dat_fcasts)

#Create threshold based off monthly clim from 1981-2010
dat_fcasts$target_month <- lubridate::month(dat_fcasts$forecast_date)
dat_fcasts$init_month <- lubridate::month(dat_fcasts$init_month)
thresh_forecasts <- dat_fcasts[dat_fcasts$forecast_date<="2010-12-31" & dat_fcasts$forecast_date>="1981-01-01",] %>% 
  group_by(ensemble, lead_month, init_month, target_month) %>% summarise_at(vars("sst_mean"),funs(mean(., na.rm=T))) 
colnames(thresh_forecasts)[5] <- "sst_clim"
ggplot(thresh_forecasts, aes(x=target_month, y=sst_clim, group=ensemble))+
  geom_point(aes(col=ensemble),shape=19)+
  facet_wrap(~init_month)+
  theme(legend.position = "bottom")

saveRDS(thresh_forecasts, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_4a/downscaled_forecast_monthly_thresholds_1980-2010.rds')








#------Step 4b: calculate DOWNSCALED FORECAST HCI-------
#Number of grid cells less than the monthly threshold

thresh_forecasts <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_4a/downscaled_forecast_monthly_thresholds_1980-2010.rds')

HCI_function_downscaledmodels <- function(x){
  fcast_area <- as.data.frame(matrix(NA, nrow=2500, ncol=5))
  colnames(fcast_area) <- c("ensemble","init_month","lead_month","forecast_date","sst_area")
  fcast_area$forecast_date <- as.Date(fcast_area$forecast_date,"%Y-%m-%d")
  fcast_area$init_month <- as.Date(fcast_area$init_month,"%Y-%m-%d")
  counter=1
  for (f in 1:length(x)){
    nc <- x[f]
    nc <- nc_open(nc)
    
    #Get vars from netcdf
    lat <- ncvar_get(nc,'lat')[1,]
    lon <- ncvar_get(nc,'lon')[,1]
    lead <- ncvar_get(nc,'lead_time')
    init_time <- as.Date(ncvar_get(nc,'init_time'), "1970-01-01")
    
    #Find dims for HCI area
    # Box: -120.3 - -116; 30.8-34.5
    lat_dimS <- which(lat==35)
    lat_dimN <- which(lat==40)
    lon_dimW <- which(lon==-127)#rough, needs to be 150km from shore
    lon_dimE <- which(lon==-120)#rough, needs to be 150km from shore
    
    #Get monthly climatology from 2003-2010
    dat <- ncvar_get(nc, 'sst')
    
    for (l in 1:12){ #lead times
      for (y in 1:29){ #years (because each forecast initialization on January OR July)
        
        #Template to force netcdf to raster dims
        template <- raster(xmn=-127,xmx=-120,ymn=35,ymx=40,resolution=0.1)
        
        #Get SST in raster form
        sst_hist <- dat[lon_dimW:lon_dimE,lat_dimS:lat_dimN,l,y]
        lat_dim <- lat[which(lat>=35 & lat<=40)]
        lon_dim <- lon[which(lon>=-127 & lon<=-120)]
        colnames(sst_hist) <- lat_dim
        rownames(sst_hist) <- lon_dim
        sst_hist <- melt(sst_hist)
        sst_r <- rasterFromXYZ(sst_hist)
        sst_r=raster::resample(sst_r,template)
        # plot(sst_r)
        # map('worldHires', add=T, fill=T)
        
        #Get distacne to shore
        d <- readRDS('~/Dropbox/PROJECTS/WRAP Location/r_distcoast.rds')
        d[d[] >= 150] = NA
        d <- crop(d,extent(-127,-120,35,40))
        # plot(d)
        # map('worldHires', add=T, fill=T)
        
        #Now crop sst to distance 
        fin <- mask(sst_r, d)
        # plot(fin)
        # map('worldHires', add=T, fill=T)
        
        
        area <- as.data.frame(rasterToPoints(fin))
        
        init <- unlist(strsplit(x[f],"_"))[5]
        init_month <- ifelse(init=="init1",1,7)
        lead.i <- lead[l] - 0.5
        forecast_date <- ymd(init_time[y]) %m+% months(lead.i)
        e <- unlist(strsplit(x[f],"_"))[4] 
        t <- thresh_forecasts$sst_clim[thresh_forecasts$ensemble==e &
                                         thresh_forecasts$lead_month == lead[l] &
                                         thresh_forecasts$init_month==init_month &
                                         thresh_forecasts$target_month==lubridate::month(forecast_date)]
        out <- nrow(area[area$value <= t,])
    
        #Write out data
        fcast_area[counter,1] <- e
        fcast_area[counter,2] <- init_time[y]
        fcast_area[counter,3] <- lead[l]  #lead month
        fcast_area[counter,4] <- forecast_date #get forecast date
        fcast_area[counter,5] <- out
        counter=counter+1
        }
      }
    }
  return(fcast_area)
}

t1 <- Sys.time()
fcast_area_init1 <- HCI_function_downscaledmodels(x=sst_files_init1)
saveRDS(fcast_area_init1,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_4b/downscaled_fcast_area_init1.rds')
fcast_area_init7 <- HCI_function_downscaledmodels(x=sst_files_init7)
saveRDS(fcast_area_init7,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_4b/downscaled_fcast_area_init7.rds')
t2 <- Sys.time()
t2-t1

summary(fcast_area_init1)
summary(fcast_area_init7)

fcast_area_init1$area <- fcast_area_init1$sst_area/max(fcast_area_init1$sst_area, na.rm=T)
fcast_area_init7$area <- fcast_area_init7$sst_area/max(fcast_area_init7$sst_area, na.rm=T)

fcast_area <- rbind(fcast_area_init1, fcast_area_init7)
fcast_area <- na.omit(fcast_area)
saveRDS(fcast_area,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_4b/downscaled_forecast_allmembers_HCI_1981-2010.rds')

fcast_area <- fcast_area[fcast_area$ensemble!="ensmean",]
fcast_area_ensemble <- fcast_area %>% group_by(init_month, lead_month, forecast_date) %>% 
  summarise_at("area", mean)
saveRDS(fcast_area_ensemble,'~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_4b/downscaled_forecast_ensmean_HCI_1981-2010.rds')


mean(fcast_area_ensemble$area, na.rm=T)
sd(fcast_area_ensemble$area, na.rm=T)
fcast_area_ensemble$init_month_month <- lubridate::month(fcast_area_ensemble$init_month)
g1 <- ggplot(data = fcast_area_ensemble[fcast_area_ensemble$init_month_month==1,], aes(x=forecast_date,y=area))+
  geom_line()+
  geom_hline(yintercept=0.37, linetype="dotted")+
  geom_hline(yintercept=0.37+0.21, col="blue")+
  geom_hline(yintercept=0.37-0.21, col="blue")+
  facet_wrap(~lead_month)+
  ggtitle("January Initialisation: downscaled forecasts")
g7 <- ggplot(data = fcast_area_ensemble[fcast_area_ensemble$init_month_month==7,], aes(x=forecast_date,y=area))+
  geom_line()+
  geom_hline(yintercept=0.37, linetype="dotted")+
  geom_hline(yintercept=0.37+0.21, col="blue")+
  geom_hline(yintercept=0.37-0.21, col="blue")+
  facet_wrap(~lead_month)+
  ggtitle("July Initialisation: downscaled forecasts")

tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_4b/init1_downscaled_forecast_HCI_timeseries_1981-2010.tiff', units="in", res=300, width=12, height=9)
plot(g1)
dev.off()

tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_4b/init7_downscaled_forecast_HCI_timeseries_1981-2010.tiff', units="in", res=300, width=12, height=9)
plot(g7)
dev.off()



#------Step 5: evaluate forecast skill of downscaled forecasts------
#Correlation; TNR; TPR: Accuracy; SEDI; Brier

#Merge in full data
fcast_area_ensemble <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_4b/downscaled_forecast_ensmean_HCI_1981-2010.rds')
obs_area <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_1b/observed_HCI_1981-2021.rds')
fcast_area_ensemble$date <- fcast_area_ensemble$forecast_date
colnames(fcast_area_ensemble)[4] <- "fcast_area"
colnames(obs_area)[3] <- "obs_area"
head(fcast_area_ensemble)
head(obs_area)
data <- left_join(fcast_area_ensemble, obs_area[,3:4], "date")
data$month <- lubridate::month(data$date)
head(data)
data$init_month_month <- lubridate::month(data$init_month)

#Plot correlations
ggplot(data[data$init_month_month==1,], aes(x=fcast_area, y=obs_area))+
  geom_point()+  
  stat_cor(method = "pearson", label.x = 0, label.y = 1, size=3)+
  facet_wrap(~lead_month)+
  ggtitle("Downscaled forecasts: January Initialisation")

ggplot(data[data$init_month_month==7,], aes(x=fcast_area, y=obs_area))+
  geom_point()+  
  stat_cor(method = "pearson", label.x = 0, label.y = 1, size=3)+
  facet_wrap(~lead_month)+
  ggtitle("Downscaled forecasts: July Initialisation")


#Compute all skill metrics
compute_metrics <- function(data=data){
  metric_summary <- as.data.frame(matrix(NA,nrow=144,ncol=9))
  colnames(metric_summary) <- c("month","lead","cor","cor_sig",
                                "tpr", "tnr","acc","sedi","bss")
  counter=1
  for (l in 0.5:11.5){ #start lead time, which then adds 6 months to get lead end time. 
    for (m in 1:12){
     
      data_temp <- data[data$lead_month==l & data$month==m,]
     
      if(nrow(data_temp)>0){
      #Correlation coefficients
      cor_i <- cor(data_temp$fcast_area, data_temp$obs_area, use="na.or.complete")
      
      #Calculate N effective degrees of freedom to determine if correlation is significant 
      N <- length(unique(lubridate::year(data_temp$forecast_date)))
      tau <- 0:(N-1)
      dm <- list()
      for (t in 1:length(tau)){
        if(t<N-1){
          rx <- cor(data_temp$fcast_area[1:(N-tau[t])],data_temp$fcast_area[(tau[t]+1):N])
          ry <- cor(data_temp$obs_area[1:(N-tau[t])],data_temp$obs_area[(tau[t]+1):N], use="na.or.complete")
        } else {
          rx <- 1
          ry <- 1
        }
        dm[[t]] <- (1-tau[t]/N) * rx * ry
      }
      dm <- unlist(dm)
      NEFF <- N/sum(dm, na.rm=T)  
      # Z <- 0.5 * log((1+cor_i)/(1-cor_i))
      # FisherZ(cor_i) #does the same as above
      # FisherZInv(Z) #just a check
      ci <- CorCI(cor_i,n=NEFF, conf.level = 0.95)
      ci_l <- ifelse(ci[2]>0,1,0 )
      } else {
        cor_i <- NA
        ci_l <- NA
      }
      
      qo = mean(data$obs_area[data$lead_month==l], na.rm=T)
      qf = mean(data$fcast_area[data$lead_month==l ], na.rm=T)
      
      #TPR: TP/TP+FN
      tpr_i <- tpr(o = data$obs_area[data$lead_month==l & data$month==m], f = data$fcast_area[data$lead_month==l & data$month==m], qo = qo, qf = qf)
      
      #TNR: TN/TN+FP
      tnr_i <- tnr(o = data$obs_area[data$lead_month==l & data$month==m], f = data$fcast_area[data$lead_month==l & data$month==m], qo = qo, qf = qf)
      
      #ACC: TN + TN / total
      acc_i <- acc(o = data$obs_area[data$lead_month==l & data$month==m], f = data$fcast_area[data$lead_month==l & data$month==m], qo = qo, qf = qf)
      
      
      # obs <- data$obs_area[data$lead_month==l & data$month==m]
      # pred <-  data$fcast_area[data$lead_month==l & data$month==m]
      # v <- verify(obs = ifelse(obs<=qo,1,0), pred = ifelse(pred<=qf,1,0), frcst.type = "binary", obs.type = "binary")
      
      
      #SEDI
      obs <- data$obs_area[data$lead_month==l & data$month==m]
      pred <-  data$fcast_area[data$lead_month==l & data$month==m]
      if (any(is.na(obs))){
        na_obs <- which(is.na(obs))
        obs <- na.omit(obs)
        pred <- pred[-na_obs]
      }
      v <- verify(obs = ifelse(obs<=qo,1,0), pred = ifelse(pred<=qf,1,0), frcst.type = "binary", obs.type = "binary")
      
      #BS and BSS
      obs_bi <- ifelse(obs <= qo,1,0)
      pred_bi <- ifelse(pred <= qf,1,0)
      BS <- mean((pred_bi - obs_bi)^2)
      BS_ref <-  0.5
      BSS <- 1 - (BS/BS_ref)
      
      metric_summary[counter,1] <- m
      metric_summary[counter,2] <- l
      metric_summary[counter,3] <- cor_i
      metric_summary[counter,4] <- ci_l
      metric_summary[counter,5] <- tpr_i
      metric_summary[counter,6] <- tnr_i
      metric_summary[counter,7] <- acc_i
      metric_summary[counter,8] <- v$SEDI
      metric_summary[counter,9] <- BSS
      counter=counter+1
    }
  }
  return(metric_summary)
}

metric_summary_init1 <- compute_metrics(data = data[data$init_month_month==1,])
metric_summary_init7 <- compute_metrics(data = data[data$init_month_month==7,])
saveRDS(metric_summary_init1, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_5/downscaled_forecast_metric_summary_init1.rds')
saveRDS(metric_summary_init7, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_5/downscaled_forecast_metric_summary_init7.rds')


#Now plot
#Correlation plots
cor <- ggplot(data = metric_summary_init1, aes(x=month,y=lead,fill=cor))+
  geom_tile(color="white")+
  geom_text(aes(label = round(cor, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.38, limit = c(0,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Correlation: July Initialisation")

#TPR plots
tpr <- ggplot(data = metric_summary_init1, aes(x=month,y=lead,fill=tpr))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tpr, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TPR") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("TPR: July Initialisation")

#TNR
tnr <- ggplot(data = metric_summary_init1, aes(x=month,y=lead,fill=tnr))+
  geom_tile(color="white")+
  geom_text(aes(label = round(tnr, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="TNR") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("TNR: July Initialisation")

#ACC
acc <- ggplot(data = metric_summary_init1, aes(x=month,y=lead,fill=acc))+
  geom_tile(color="white")+
  geom_text(aes(label = round(acc, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0.5, limit = c(0,1), space = "Lab",
                       name="Accuracy") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Accuracy: July Initialisation")
# (0.21*0.21) + (0.79*0.79) = 0.67
# (0.5*0.5) + (0.5*0.5) = 0.5 #based on median (50% quantile) being roughly 0.5

#SEDI
sedi <- ggplot(data = metric_summary_init1, aes(x=month,y=lead,fill=sedi))+
  geom_tile(color="white")+
  geom_text(aes(label = round(sedi, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="SEDI") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("SEDI: July Initialisation")

bss <- ggplot(data = metric_summary_init1, aes(x=month,y=lead,fill=bss))+
  geom_tile(color="white")+
  geom_text(aes(label = round(bss, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(0,1), space = "Lab",
                       name="BSS") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("BSS: July Initialisation")

tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_5/fig_skill_plots_january_initialisation.tiff', res=200, units="in",width=15,height=8)
grid.arrange(cor, tpr, tnr, acc, sedi,bss, nrow=2)
dev.off()

#-----Step 6: compare downscaled skill with global skill------
dsl1 <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_5/downscaled_forecast_metric_summary_init1.rds')
dsl7 <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_5/downscaled_forecast_metric_summary_init7.rds')
gbl <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_3b/global_forecasts_metric_summary.rds')

dsl1 <- na.omit(dsl1)
dsl7 <- na.omit(dsl7)
dsl <- rbind(dsl1, dsl7)
dsl$model <- "downscaled"
gbl$model <- "global"
all <- rbind(dsl, gbl)

cor_long <- all %>% pivot_wider(id_cols=1:2,names_from = model,values_from=cor)
tpr_long <- all %>% pivot_wider(id_cols=1:2,names_from = model,values_from=tpr)
tnr_long <- all %>% pivot_wider(id_cols=1:2,names_from = model,values_from=tnr)
acc_long <- all %>% pivot_wider(id_cols=1:2,names_from = model,values_from=acc)
sedi_long <- all %>% pivot_wider(id_cols=1:2,names_from = model,values_from=sedi)
bss_long <- all %>% pivot_wider(id_cols=1:2,names_from = model,values_from=bss)

cor_long$diff <- cor_long$global-cor_long$downscaled
tpr_long$diff <- tpr_long$global-tpr_long$downscaled
tnr_long$diff <- tnr_long$global-tnr_long$downscaled
acc_long$diff <- acc_long$global-acc_long$downscaled
sedi_long$diff <- sedi_long$global-sedi_long$downscaled
bss_long$diff <- bss_long$global-bss_long$downscaled

g1 <- ggplot(cor_long, aes(x=month, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.5,0.5), space = "Lab", 
                       name="") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Differece in Correlation: global - downscaled")

g2 <- ggplot(tpr_long, aes(x=month, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.3,0.3), space = "Lab",
                       name="TPR") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Differece in TNR: global - downscaled")

g3 <- ggplot(tnr_long, aes(x=month, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.4,0.4), space = "Lab",
                       name="TNR") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Differece in TPR: global - downscaled")

g4 <- ggplot(acc_long, aes(x=month, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.3,0.3), space = "Lab",
                       name="Accuracy") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Differece in ACC: global - downscaled")

g5 <- ggplot(sedi_long, aes(x=month, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.8,0.8), space = "Lab",
                       name="SEDI") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Differece in SEDI: global - downscaled")

g6 <- ggplot(bss_long, aes(x=month, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 1))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-0.8,0.8), space = "Lab",
                       name="BSS") +
  scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Differece in BSS: global - downscaled")


tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_6/fig_skill_plots_comparison.tiff', res=200, units="in",width=15,height=8)
grid.arrange(g1, g2, g3, g4, g5,g6, nrow=2)
dev.off()

#-----Step 9: early warning system example of novel years-----
#Merge in full data: global model all years
fcast_area_ensemble <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/forecast_ensmean_HCI_1981-2021.rds')
obs_area <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_1b/observed_HCI_1981-2021.rds')
fcast_area_ensemble$date <- fcast_area_ensemble$forecast_date
colnames(fcast_area_ensemble)[4] <- "fcast_area"
colnames(obs_area)[3] <- "obs_area"
head(fcast_area_ensemble)
head(obs_area)
data <- left_join(fcast_area_ensemble, obs_area[,3:4], "date")
data$month <- lubridate::month(data$date)
head(data)

#When exactly did HCI occurr in 2005, and 2015
plot(data$date,data$obs_area, type='b')
ggplot(data, aes(x=date,y=obs_area))+
  geom_line()+
  geom_hline(yintercept = 0.4, linetype="dotted")+
  geom_hline(yintercept = 0.4-0.28, col="blue")+
  geom_hline(yintercept = 0.4+0.28, col="blue")
data_high <- data[data$obs_area<=0.12,]
data_high <- data_high[data_high$date>"2000-01-01",]
unique(data_high$date)

#Observed high HCI:
# [1] "2003-01-01" "2003-08-01" "2005-03-01" "2005-05-01" "2006-06-01" "2009-06-01" "2012-11-01"
# [8] "2014-07-01" "2014-09-01" "2014-10-01" "2014-11-01" "2014-12-01" "2015-01-01" "2015-02-01"
# [15] "2015-03-01" "2015-07-01" "2015-08-01" "2015-10-01" "2015-11-01" "2015-12-01" "2016-01-01"
# [22] "2016-02-01" "2016-03-01" "2016-04-01" "2016-05-01" "2016-11-01" "2017-09-01" "2017-11-01"
# [29] "2018-01-01" "2018-12-01" "2020-05-01" NA         

#Moving ahead with two examples: March & May 2005; and Sep 2014 - May 2016
thresholds <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2c/global_forecast_hci_thresholds.rds')
thresholds$high_hci <- round(thresholds$high_hci,2)
cs1 <- data
cs1 <- left_join(cs1, thresholds[,c(1,4)], by="lead_month")
cs1_march <- cs1[cs1$forecast_date=="2005-03-01",]
cs1_may <- cs1[cs1$forecast_date=="2005-05-01",]
cs1_march$hit <- ifelse(cs1_march$fcast_area<=cs1_march$high_hci,1,0)
cs1_may$hit <- ifelse(cs1_may$fcast_area<=cs1_may$high_hci,1,0)

cs2 <- data
cs2 <- left_join(cs2, thresholds[,c(1,4)], by="lead_month")
cs2 <- cs2[cs2$forecast_date>="2014-09-01" & cs2$forecast_date<="2016-05-01",]
cs2$hit <- ifelse(cs2$fcast_area<=cs2$high_hci,1,0)

cs2[cs2$date=="2014-09-01",]
view(cs2[cs2$lead_month==11.5,])

