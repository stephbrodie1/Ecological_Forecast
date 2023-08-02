#Code used to generate figures
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
library(reshape2)
library(cmocean)
library(scales)

#-----Figure 1: Maps of forecasts-------
#make maps of global and downscaled forecasts
#Decided to pick CanCM4 ens 2 for 2010-01-01

#Global model
nc <- nc_open('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/NMME_CCS_SST_forecasts/sst_forecast_CanCM4i.nc')# d="1981-01-01" #10 ens

#Get vars from netcdf
lat <- ncvar_get(nc,'lat')
lon <- ncvar_get(nc,'lon')
lead <- ncvar_get(nc,'lead')
time <- ncvar_get(nc,'time') #months since 1960-01-01. I don't know how to convert this, so next few lines are a workaround
t1 <- ymd("1981-01-01")
time <- time-time[1]
init_time <- t1 %m+% months(time)
clim_idx <- which(init_time>="1991-01-01" & init_time<="2020-12-31")
clim_start <- clim_idx[1]
clim_end <- tail(clim_idx,n=1)
member <- ncvar_get(nc,'member')
t_start <- which(init_time>="2010-01-01")[1]
lat_dimS <- which(lat==30)
lat_dimN <- which(lat==48)
lon_dimW <- which(lon==249) 
lon_dimE <- which(lon==229) 
dat <- ncvar_get(nc, 'sst')

#Manipulate data
d <- dat[lon_dimW:lon_dimE,lat_dimS:lat_dimN,2,1,t_start]
lat_dim <- lat[which(lat>=30 & lat<=48)]
lon_dim <- lon
colnames(d) <- lat_dim
rownames(d) <- lon_dim
d <- melt(d)
dr <- rasterFromXYZ(d)
dr2 <- flip(dr,1) 
plot(dr2)
dr3 <- as.data.frame(rasterToPoints(dr2))
dr3$lon <- dr3$x - 360


#get dist-to-shore layers
# d <- readRDS('~/Dropbox/PROJECTS/Completed_Projects/WRAP Location/r_distcoast.rds')
# ex <- c(-130.5,-115,35,40)
# d <- crop(d,ex)
# t75 <- rasterToContour(d,levels=c(75))
# t150 <- rasterToContour(d,levels=c(150))
# t_f75 <- fortify(t75)
# t_f150 <- fortify(t150)
# saveRDS(t_f75, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/contour_75km.rds')
# saveRDS(t_f150, '~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/contour_150km.rds')
# plot(d)
# plot(t75, add=T)
# plot(t150, add=T)

#loggerhead closure area
lc <- readShapePoly("~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/reshapefileforloggerheadclosurearea/LoggerheadClosure.shp")
lc.df <- fortify(lc)

#Plot
t_f75 <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/contour_75km.rds')
t_f150 <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/contour_150km.rds')
g1 <- ggplot()+
  geom_tile(data=dr3, aes(x=lon, y=y, fill=value))+
  geom_path(data=t_f75, aes(x=long, y=lat, group=group),col='blue', linetype="solid") + #HCI
  geom_path(data=t_f150, aes(x=long, y=lat, group=group),col="green", linetype="solid") + #HCI
  geom_segment(aes(x = -121.6837, y = 35.05, xend = -119, yend = 35.05, colour = "segment"), col="blue")+ #HCI #-121.6837 & 35.05000
  geom_segment(aes(x = -125.1444, y = 39.95, xend = -119, yend = 39.95, colour = "segment"), col="blue")+ #HCI-125.1444 &  39.95000
  geom_segment(aes(x = -122.7562, y = 35.05, xend = -121.6837, yend = 35.05, colour = "segment"), col="green")+ #HCI -122.7562 &  35.05000
  geom_segment(aes(x = -126.0797, y = 39.95, xend = -125.1444, yend = 39.95, colour = "segment"), col="green")+ #HCI -126.0797 &  39.95000
  theme_classic() +  labs(y="", x="", title = "") + 
  theme(legend.position="right",legend.title = element_blank(),rect = element_rect(fill = "transparent"))+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("thermal")(256), name="Sea Surface Temperature",limits=c(7.5,19)) +
  geom_segment(aes(x = -120, y = 31, xend = -115, yend = 31, colour = "segment"), col="black")+ #TOTAL
  geom_segment(aes(x = -120, y = 31, xend = -120, yend = 35, colour = "segment"), col="black")+ #TOTAL
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-127.5,-115),ylim=c(30.5,41.5)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_polygon(data=lc.df, aes(long, lat, group = group),color="white", fill=NA, linetype='dashed')
tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/Fig1_map_global.tiff', res=300, units="in", width=5, height=4,bg = "transparent")
plot(g1)
dev.off()

#Downscaled model
nc <- nc_open('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/ForecastFields/sst_forecast_ens2_init1_monthavg_1982-2010.nc') 
#Get vars from netcdf
lat <- ncvar_get(nc,'lat')[1,]
lon <- ncvar_get(nc,'lon')[,1]
lead <- ncvar_get(nc,'lead_time')
init_time <- as.Date(ncvar_get(nc,'init_time'), "1970-01-01")
lat_dimS <- which(lat==30)
lat_dimN <- which(lat==48)
lon_dimW <- which(lon==-134)#rough, needs to be 75km from shore
lon_dimE <- which(lon==-115.5)#rough, needs to be 75km from shore
dat <- ncvar_get(nc, 'sst')
sst_hist <- dat[lon_dimW:lon_dimE,lat_dimS:lat_dimN,1,29]
lat_dim <- lat
lon_dim <- lon
colnames(sst_hist) <- lat_dim
rownames(sst_hist) <- lon_dim
sst_hist <- melt(sst_hist)
sst_r <- rasterFromXYZ(sst_hist)
sst_r <- as.data.frame(rasterToPoints(sst_r))
g2 <- ggplot()+
  geom_tile(data=sst_r, aes(x=x, y=y, fill=value))+
  geom_path(data=t_f75, aes(x=long, y=lat, group=group),col='blue', linetype="solid") + #HCI
  geom_path(data=t_f150, aes(x=long, y=lat, group=group),col="green", linetype="solid") + #HCI
  geom_segment(aes(x = -121.6837, y = 35.05, xend = -119, yend = 35.05, colour = "segment"), col="blue")+ #HCI #-121.6837 & 35.05000
  geom_segment(aes(x = -125.1444, y = 39.95, xend = -119, yend = 39.95, colour = "segment"), col="blue")+ #HCI-125.1444 &  39.95000
  geom_segment(aes(x = -122.7562, y = 35.05, xend = -121.6837, yend = 35.05, colour = "segment"), col="green")+ #HCI -122.7562 &  35.05000
  geom_segment(aes(x = -126.0797, y = 39.95, xend = -125.1444, yend = 39.95, colour = "segment"), col="green")+ #HCI -126.0797 &  39.95000
  theme_classic() +  labs(y="", x="", title = "") + 
  theme(legend.position="right",legend.title = element_blank(),rect = element_rect(fill = "transparent"))+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("thermal")(256), name="Sea Surface Temperature",limits=c(7.5,19)) +
  geom_segment(aes(x = -120, y = 31, xend = -115, yend = 31, colour = "segment"), col="black")+
  geom_segment(aes(x = -120, y = 31, xend = -120, yend = 35, colour = "segment"), col="black")+
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-127.5,-115),ylim=c(30.5,41.5)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_polygon(data=lc.df, aes(long, lat, group = group),color="white", fill=NA, linetype='dashed')
tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/Fig1_map_downscaled.tiff', res=300, units="in", width=5, height=4, bg = "transparent")
plot(g2)
dev.off()
    
grid.arrange(g1, g2, nrow=1)


#-----Figure 1: timeseries of tools-----
#HCI:
hc <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_1b/observed_HCI_1981-2021.rds')
g1 <- ggplot(hc, aes(x=date,y=area))+
  geom_line()+
  # geom_segment(aes(x = as.Date("1981-01-01"), y = 0.4, xend = as.Date("2020-12-31"), yend = 0.4),linetype="dashed")+ 
  # annotate( "rect",xmin=as.Date("1981-01-01"), xmax=as.Date("2020-12-31"), alpha=0.5, ymin=0, ymax=0.4, fill="red")+
  # annotate( "rect",xmin=as.Date("1981-01-01"), xmax=as.Date("2020-12-31"), alpha=0.5, ymin=0.4+0.28, ymax=Inf, fill="blue")+
  # annotate( "rect",xmin=as.Date("1981-01-01"), xmax=as.Date("2020-12-31"), alpha=0.5, ymin=0.4, ymax=0.4+0.28, fill="orange")+
  geom_hline(yintercept = 0.4,col="red")+
  scale_x_date(expand=c(0,0),limits = c(as.Date("1981-01-01"),as.Date("2020-12-31")),breaks = scales::breaks_pretty(10)) +
  theme_classic()+
  ylab("HCI") + xlab("") +
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1))+  #makes a box
  theme(legend.position = "none",rect = element_rect(fill = "transparent"))
  # scale_x_continuous(breaks = round(seq(as.Date("1981-01-01"), as.Date("1981-01-01"), by = 365),1))
tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/Fig1_timeseries_HCI.tiff', units="in", res=300, width=6, height=2, bg = "transparent")
plot(g1)
dev.off()


tt <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1a/observed_anomalies_bestcase.rds')
readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1a/observed_threshold_bestcase.rds')

g2 <- ggplot(tt, aes(x=date,y=roll_obs_anom))+
  geom_line()+
  geom_hline(yintercept=0.41, col="red")+
  scale_x_date(expand=c(0,0),limits = c(as.Date("1981-01-01"),as.Date("2020-12-31")),breaks = scales::breaks_pretty(10)) +
  theme_classic()+
  ylab("SSTA") + xlab("") +
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1))+  #makes a box
  theme(legend.position = "none",rect = element_rect(fill = "transparent"))
tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/Fig1_timeseries_TOTAL.tiff', units="in", res=300, width=6, height=2, bg = "transparent")
plot(g2)
dev.off()


#-----Figure 2: HCI time series and skill assessment------
#HCI time series:
#HCI time series:
hc_obs <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_1b/observed_HCI_1981-2021.rds')
hc_fc <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/forecast_ensmean_HCI_1981-2021.rds')
th <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2c/global_forecast_hci_thresholds.rds')

o <- hc_obs
o$forecast_date <- o$date
o$type = "Observed"
o <- o[,c(5,6,3)]
o$lead_month <- 0
o$mean <- 0.4 #threshold

fm <- left_join(hc_fc,th[,1:2],"lead_month")
f <- fm[fm$lead_month==0.5 |  fm$lead_month==6.5 | fm$lead_month==11.5,]
# f <- fm
f$type = "Forecast"
f <- f[,c(3,6,4,2,5)]
dt <- rbind(o,f)

dt$type <- as.factor(dt$type)
dt$type <- fct_rev(dt$type)

dt <- dt[dt$forecast_date>="2013-01-01",]
dt <- dt[dt$forecast_date<="2017-12-31",]
dt$lead_month <- as.factor(dt$lead_month)
dt <- na.omit(dt)

dt$group <- ifelse(dt$area<=dt$mean,1,0)

g1 <- ggplot()+
  geom_line(data= dt, aes(x=forecast_date,y=area, color=as.factor(lead_month)),size=0.7)+
  geom_hline(data= dt[dt$type=="Forecast",], aes(yintercept=mean, color=lead_month))+
  geom_hline(yintercept=0.4, color="black")+
  geom_rect(aes(xmin=as.Date("2014-03-01"), xmax=as.Date("2016-12-31"), ymin=-0.15, ymax=Inf),alpha=0.1)+
  scale_x_date(expand=c(0,0),limits = c(as.Date("2013-01-01"),as.Date("2017-12-31")),breaks = scales::breaks_pretty(10)) +
  scale_y_continuous(expand=c(0,0),limits=c(-0.15,1.15)) +
  theme_classic()+
  ylab("HCI") + xlab("") + 
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1))+  #makes a box
  theme(legend.position = "bottom", legend.title = element_blank())+
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Habitat Compression Index")+
  scale_color_manual(values = c("black","#F8766D", "#00BA38" ,"#619CFF"),labels = c("Observed","0.5", "6.5", "11.5"))


#poster
# g_poster <- ggplot()+
#   geom_line(data= dt, aes(x=forecast_date,y=area, color=as.factor(lead_month)),size=0.7)+
#   geom_hline(data= dt[dt$type=="Forecast",], aes(yintercept=mean, color=lead_month))+
#   geom_hline(yintercept=0.4, color="black")+
#   geom_rect(aes(xmin=as.Date("2014-03-01"), xmax=as.Date("2016-12-31"), ymin=-0.15, ymax=Inf),alpha=0.2)+
#   scale_x_date(expand=c(0,0),limits = c(as.Date("2013-01-01"),as.Date("2017-12-31")),breaks = scales::breaks_pretty(10)) +
#   scale_y_continuous(expand=c(0,0),limits=c(-0.15,1.15)) +
#   theme_classic()+
#   ylab("HCI") + xlab("") + 
#   theme( panel.border = element_rect(colour = "black", fill=NA, size=1))+  #makes a box
#   theme(legend.position = "bottom", legend.title = element_blank())+
#   theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   theme(rect = element_rect(fill = "transparent"))+
#   ggtitle("Habitat Compression Index")+
#   scale_color_manual(values = c("black","#F8766D", "#00BA38" ,"#619CFF"),labels = c("Observed","0.5", "6.5", "11.5"))
# tiff('~/Dropbox/Conference and Workshops/AMSA 2023/hci_fig.tiff', res=300, height=4, width = 7, units="in", bg="transparent")
# plot(g_poster)
# dev.off()


# OLD TIMESERIES
# hc_obs <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_1b/observed_HCI_1981-2021.rds')
# hc_fc <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2b/forecast_ensmean_HCI_1981-2021.rds')
# th <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_2c/global_forecast_hci_thresholds.rds')
# 
# o <- hc_obs
# o$forecast_date <- o$date
# o$type = "Observed"
# o <- o[,c(5,6,3)]
# f <- hc_fc[hc_fc$lead_month==0.5,]
# f$type = "Forecast"
# f <- f[,c(3,5,4)]
# dt <- rbind(o,f)
# 
# dt$type <- as.factor(dt$type)
# dt$type <- fct_rev(dt$type)
# g1 <- ggplot(dt, aes(x=forecast_date,y=area, color=type))+
#   geom_line()+
#   # scale_color_manual(values=c('#e74c3c','#2c3e50'))+
#   annotate( "rect",xmin=as.Date("1981-01-01"), xmax=as.Date("2020-12-31"), alpha=0.2, ymin=0, ymax=0.419, fill="#00BFC4")+#, col="#00BFC4")+
#   annotate( "rect",xmin=as.Date("1981-01-01"), xmax=as.Date("2020-12-31"), alpha=0.2, ymin=0, ymax=0.4, fill="#F8766D")+#, col="#F8766D")+
#   # annotate( "rect",xmin=as.Date("1981-01-01"), xmax=as.Date("2020-12-31"), alpha=0.2, ymin=0.353, ymax=0.445, fill="red", col="red")+
#   scale_x_date(expand=c(0,0),limits = c(as.Date("1981-01-01"),as.Date("2020-12-31")),breaks = scales::breaks_pretty(10)) +
#   theme_classic()+
#   ylab("HCI") + xlab("") +
#   theme( panel.border = element_rect(colour = "black", fill=NA, size=1))+  #makes a box
#   theme(legend.position = "bottom", legend.title = element_blank())

##skill maps
metric_summary <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_3a/global_forecasts_metric_summary.rds')
# persistence <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_7b/persistence_forecast_metric_summary.rds')
# persistence$lead <- persistence$lead-0.5

# cors <- left_join(metric_summary[,c(1,2,3,4)],persistence[,c(1,2,3,4)],by=c("month","lead"))
# cors$diff <- cors$cor.x - cors$cor.y
# 
# accs <- left_join(metric_summary[,c(1,2,7)],persistence[,c(1,2,7)],by=c("month","lead"))
# accs$diff <- accs$acc.x - accs$acc.y

#Now plot
#Correlation plots
cor <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=cor))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = cmocean("balance")(256),
                        values=rescale(c(-1,0,1)),limits=c(-1,1),name="")+
  geom_point(data = metric_summary[metric_summary$cor_sig==1,], col="black", size=0.5)+
  scale_x_continuous(breaks = seq_along(month.abb), labels = month.abb, expand = c(0, 0))+
  scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
  labs(x="Forecast Month", y="Lead time (months)")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.5,'cm'))+
  ggtitle("Correlation")
#ACC
acc <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=acc))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(0,0.5,1)),limits=c(0,1),name="")+
  scale_x_continuous(breaks = seq_along(month.abb), labels = month.abb, expand = c(0, 0))+
  scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
  labs(x="Forecast Month", y="Lead time (months)")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.5,'cm'))+
  ggtitle("Accuracy")
# (0.21*0.21) + (0.79*0.79) = 0.67
# (0.5*0.5) + (0.5*0.5) = 0.5 #based on median (50% quantile) being roughly 0.5

#SEDI
sedi <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=sedi))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(-0.5,0,1)),limits=c(-0.5,1),name="")+
  scale_x_continuous(breaks = seq_along(month.abb), labels = month.abb, expand = c(0, 0))+
  scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
  labs(x="Forecast Month", y="Lead time (months)")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(0.5,'cm'))+
  ggtitle("SEDI")


tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/Fig2_HCI_timeseries_and_skill.tiff', units="in", res=300, width=9, height=7)
ggarrange(g1,  # First row with scatter plot
          ggarrange(cor, acc, sedi, ncol = 3, labels = c("B", "C", "D")), # Second row with box and dot plots
          nrow = 2, labels = "A" ) 
dev.off()


#POSTER plot

#Now plot
# #Correlation plots
# cor_poster <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=cor))+
#   geom_tile(color="white")+
#   scale_fill_gradientn(colors = cmocean("balance")(256),
#                        values=rescale(c(-1,0,1)),limits=c(-1,1),name="")+
#   geom_point(data = metric_summary[metric_summary$cor_sig==1,], col="black", size=0.5)+
#   scale_x_continuous(breaks = seq_along(month.abb), labels = month.abb, expand = c(0, 0))+
#   scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
#   labs(x="Forecast Month", y="Lead time (months)")+
#   theme(legend.position = "bottom",
#         axis.text.x = element_text(angle = 45),
#         legend.key.width = unit(1, "cm"),
#         legend.key.height = unit(0.5,'cm'))+
#   theme(rect = element_rect(fill = "transparent"))+
#   ggtitle("Correlation") 
# tiff('~/Dropbox/Conference and Workshops/AMSA 2023/cor_fig.tiff', res=300, height=4, width = 4, units="in", bg="transparent")
# plot(cor_poster)
# dev.off()



#-----Figure 3: TOTAL time series and skill assessment-----

#TOTAL timeseries
tt_obs <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1a/observed_anomalies_bestcase.rds')
th_obs <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1a/observed_threshold_bestcase.rds')
tt_fc <- readRDS("~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2a/global_fcasts.rds")
th_fc <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2a/thresholds_globalforecasts_bestcase.rds')

#get 6-month rolling mean of forecast
ens <- tt_fc %>% group_by(lead_month, init_month, forecast_date) %>% summarise_at("sst_anom",mean, na.rm=T)
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
global_ensemble_rolling$month <- lubridate::month(global_ensemble_rolling$init_month)

#Merge observed and forecasts
global_ensemble_rolling$forecast_date <- global_ensemble_rolling$init_month# %m+%  months(global_ensemble_rolling$lead+0.5)
global_ensemble_rolling$type = "Forecast"
global_ensemble_rolling <- left_join(global_ensemble_rolling,th_fc,by= "lead")
tt_obs$type = "Observed"
tt_obs$forecast_date = tt_obs$date
tt_obs$sst_anom_rolling <- tt_obs$roll_obs_anom
tt_obs$ssta_quantile <- 0.41

f <- global_ensemble_rolling[,c(5,6,3,2,7)]
tt_obs$lead <- 0 
o <- tt_obs[,c(10,9,11,13,12)]
dt <- rbind(o,f)

dt$lead <- dt$lead+6
dt$lead <- as.factor(dt$lead)

dt <- dt[dt$forecast_date>="2013-01-01",]
dt <- dt[dt$forecast_date<="2017-12-31",]
dt <- na.omit(dt)

dt <- dt[dt$lead==6 |dt$lead==6.5 | dt$lead==8.5 | dt$lead==11.5 ,]

g1 <- ggplot()+
  geom_line(data = dt, aes(x=forecast_date,y=sst_anom_rolling, color=as.factor(lead)),size=0.7)+
  geom_hline(data = dt[dt$type=="Forecast",], aes(yintercept = ssta_quantile, color=lead))+
  geom_hline(yintercept = 0.41, color="black")+
  geom_rect(aes(xmin=as.Date("2014-03-01"), xmax=as.Date("2016-12-31"), ymin=-0.4, ymax=Inf),alpha=0.1)+
  scale_x_date(expand=c(0,0),limits = c(as.Date("2013-01-01"),as.Date("2017-12-31")),breaks = scales::breaks_pretty(10)) +
  scale_y_continuous(expand = c(0, 0),limits=c(-0.4,2)) +
  theme_classic()+
  ylab("SSTA") + xlab("") +
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1))+  #makes a box
  theme(legend.position = "bottom", legend.title = element_blank())+
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Temperature Observations to Avoid Loggerheads")+
  scale_color_manual(values = c("black","#F8766D", "#00BA38" ,"#619CFF"),labels = c("Observed","6.5", "8.5","11.5"))


#OLD time series code
# #TOTAL time series:
# tt_obs <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1a/observed_anomalies_bestcase.rds')
# th_obs <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_1a/observed_threshold_bestcase.rds')
# tt_fc <- readRDS("~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2a/global_fcasts.rds")
# th_fc <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_2a/thresholds_globalforecasts_bestcase.rds')
# 
# #get 6-month rolling mean of forecast
# ens <- tt_fc %>% group_by(lead_month, init_month, forecast_date) %>% summarise_at("sst_anom",mean, na.rm=T)
# global_ensemble_rolling <- as.data.frame(matrix(NA, nrow=100,ncol=3))
# colnames(global_ensemble_rolling) <- c("init_month","lead","sst_anom_rolling")
# global_ensemble_rolling$init_month <- as.Date(global_ensemble_rolling$init_month)
# counter=1
# dates <- as.Date(unique(ens$init_month),"%Y-%m-$d")
# for (d in 1:length(dates)){
#   for (l in 0.5:5.5){
#     dt <- dates[d]
#     sst_anom_fcast_rolling <- mean(ens$sst_anom[ens$init_month==dt &
#                                                   ens$lead_month>=l & ens$lead_month<=l+6],na.rm=T)
#     global_ensemble_rolling[counter,1] <- as.Date(dt)
#     global_ensemble_rolling[counter,2] <- l
#     global_ensemble_rolling[counter,3] <- sst_anom_fcast_rolling
#     counter=counter+1
#   }
# }
# global_ensemble_rolling$month <- lubridate::month(global_ensemble_rolling$init_month)
# 
# 
# #Merge observed and forecasts
# global_ensemble_rolling$forecast_date <- global_ensemble_rolling$init_month #%m+%  months(global_ensemble_rolling$lead+0.5)
# global_ensemble_rolling <- global_ensemble_rolling[global_ensemble_rolling$lead==0.5,] 
# 
# global_ensemble_rolling$type = "Forecast"
# tt_obs$type = "Observed"
# tt_obs$forecast_date = tt_obs$date
# tt_obs$sst_anom_rolling <- tt_obs$roll_obs_anom
# f <- global_ensemble_rolling[,c(5,6,3)]
# o <- tt_obs[,c(10,9,11)]
# dt <- rbind(o,f)
# 
# dt$type <- as.factor(dt$type)
# dt$type <- fct_rev(dt$type)
# g1 <- ggplot(dt, aes(x=forecast_date,y=sst_anom_rolling, color=type))+
#   geom_line()+
#   geom_hline(yintercept=0.12,col="#00BFC4")+
#   geom_hline(yintercept=0.41,col="#F8766D")+
#   scale_x_date(expand=c(0,0),limits = c(as.Date("1981-01-01"),as.Date("2020-12-31")),breaks = scales::breaks_pretty(10)) +
#   theme_classic()+
#   ylab("SSTA") + xlab("") +
#   theme( panel.border = element_rect(colour = "black", fill=NA, size=1))+  #makes a box
#   theme(legend.position = "bottom", legend.title = element_blank())


##skill maps
metric_summary <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_8a/metric_summary_global_1980-2020_bestcasescenario.rds')

# cors <- left_join(metric_summary[,c(1,2,3)],persistence[,c(1,2,3)],by=c("month","lead"))
# cors$diff <- cors$cor.x - cors$cor.y
# accs <- left_join(metric_summary[,c(1,2,6)],persistence[,c(1,2,6)],by=c("month","lead"))
# accs$diff <- accs$acc.x - accs$acc.y

#Now plot
#Correlation plots
metric_summary$correlation_months <- as.factor(metric_summary$correlation_months)
metric_summary$correlation_months <- factor(metric_summary$correlation_months, levels=rev(levels(metric_summary$correlation_months)))
metric_summary$lead <- metric_summary$lead + 6 
cor <- ggplot(data = metric_summary, aes(x=correlation_months,y=lead,fill=correlation_values))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(-1,0,1)),limits=c(-1,1),name="")+
  geom_point(data = metric_summary[metric_summary$sig_values==1,], col="black", size=0.5)+
  scale_x_discrete(expand = c(0, 0), labels=c("Jun","Jul","Aug")) +
  scale_y_continuous(breaks=6.5:11.5, expand = c(0, 0)) +
  labs(x="Closure Month", y="Lead time (months)")+
  theme(legend.position = "bottom",
        legend.key.width = unit(0.8, "cm"),
        legend.key.height = unit(0.4,'cm'))+
  ggtitle("Correlation")

#ACC
metric_summary$acc_months <- as.factor(metric_summary$acc_months)
metric_summary$acc_months <- factor(metric_summary$acc_months, levels=rev(levels(metric_summary$acc_months)))
acc <- ggplot(data = metric_summary, aes(x=acc_months,y=lead,fill=acc_values))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(0,0.62,1)),limits=c(0,1),name="")+
  scale_x_discrete(expand = c(0, 0), labels=c("Jun","Jul","Aug")) +
  scale_y_continuous(breaks=6.5:11.5, expand = c(0, 0)) +
  labs(x="Closure Month", y="Lead time (months)")+
  theme(legend.position = "bottom",
        legend.key.width = unit(0.8, "cm"),
        legend.key.height = unit(0.4,'cm'))+
  ggtitle("Accuracy")

#SEDI
metric_summary$sedi_months <- as.factor(metric_summary$sedi_months)
metric_summary$sedi_months <- factor(metric_summary$sedi_months, levels=rev(levels(metric_summary$sedi_months)))
sedi <- ggplot(data = metric_summary, aes(x=sedi_months,y=lead,fill=sedi_values))+
  geom_tile(color="white")+
  # geom_text(aes(label = round(sedi, 1))) +
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(-0.5,0,1)),limits=c(-0.5,1),name="")+
  scale_x_discrete(expand = c(0, 0), labels=c("Jun","Jul","Aug")) +
  scale_y_continuous(breaks=6.5:11.5, expand = c(0, 0)) +
  labs(x="Closure Month", y="Lead time (months)")+
  theme(legend.position = "bottom",
        legend.key.width = unit(0.8, "cm"),
        legend.key.height = unit(0.4,'cm'))+
  ggtitle("SEDI")


tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/Fig3_TOTAL_timeseries_and_skill.tiff', units="in", res=300, width=7, height=6)
ggarrange(g1,  # First row with scatter plot
          ggarrange(cor, acc, sedi, ncol = 3, labels = c("B", "C", "D")), # Second row with box and dot plots
          nrow = 2, labels = "A" ) 
dev.off()

#-----Figure 4: Boxplots of scenario performance-----
#Boxplot of HCI and TOTAL, each with 3 boxes (downscaled skill, global reduced skill, and global comparative skill)

#HCI
metric_summary <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_3b/global_forecasts_metric_summary.rds')
metric_summary_init1 <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_5/downscaled_forecast_metric_summary_init1.rds')
metric_summary_init7 <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_5/downscaled_forecast_metric_summary_init7.rds')
metric_summary_init1 <- na.omit(metric_summary_init1)
metric_summary_init7 <- na.omit(metric_summary_init7)
ds <- rbind(metric_summary_init1, metric_summary_init7)
metric_summary_r <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_3c/global_forecasts_metric_summary_only3.rds')

#acc
d <- left_join(ds[,c(1,2,6)], metric_summary[,c(1,2,6)], c("month","lead"))
colnames(d)[3:4] <- c("Downscaled","Global Full")
d <- left_join(d, metric_summary_r[,c(1,2,6)], c("month","lead"))
colnames(d)[5] <- c("Global Reduced")
data_d <- d %>% pivot_longer(c(3:5),names_to = "model", values_to = "values")
acc <- data_d
acc$tool <- "HCI"
boxplot(values~model,data=acc)

#cor
d <- left_join(ds[,c(1,2,3)], metric_summary[,c(1,2,3)], c("month","lead"))
colnames(d)[3:4] <- c("Downscaled","Global Full")
d <- left_join(d, metric_summary_r[,c(1,2,3)], c("month","lead"))
colnames(d)[5] <- c("Global Reduced")
data_d <- d %>% pivot_longer(c(3:5),names_to = "model", values_to = "values")
cor <- data_d
cor$tool <- "HCI"
boxplot(values~model,data=cor)

#sedi
d <- left_join(ds[,c(1,2,7)], metric_summary[,c(1,2,7)], c("month","lead"))
colnames(d)[3:4] <- c("Downscaled","Global Full")
d <- left_join(d, metric_summary_r[,c(1,2,7)], c("month","lead"))
colnames(d)[5] <- c("Global Reduced")
data_d <- d %>% pivot_longer(c(3:5),names_to = "model", values_to = "values")
sedi <- data_d
sedi$tool <- "HCI"
boxplot(values~model,data=acc)




#TOTAL
##skill maps for comparative scenario (i.e. 73 ensemble members)
gl <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_8b/metric_summary_global_1980-2010_applesscenario.rds')
ds <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_9b/metric_downscaled.rds')
ds <- na.omit(ds)
gl_r <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_8c/metric_summary_global_1980-2010_only3_scenario.rds')


#acc
b <- left_join(ds[,c(1,8,9)], gl[,c(1,8,9)], c("lead", "acc_months"))
colnames(b)[3:4] <- c("Downscaled", "Global Full")
b <- left_join(b, gl_r[,c(1,8,9)], c("lead", "acc_months"))
colnames(b)[5] <- c("Global Reduced")
colnames(b)[2] <- c("month")
b <- b[!duplicated(b),]
data_b <- b %>% pivot_longer(c(3:5),names_to = "model", values_to = "values")
acc_total <- data_b
acc_total$tool <- "TOTAL"
boxplot(values~model,data=acc_total)

#cor
b <- left_join(ds[,c(1,2,3)], gl[,c(1,2,3)], c("lead", "correlation_months"))
colnames(b)[3:4] <- c("Downscaled", "Global Full")
b <- left_join(b, gl_r[,c(1,2,3)], c("lead", "correlation_months"))
colnames(b)[5] <- c("Global Reduced")
colnames(b)[2] <- c("month")
b <- b[!duplicated(b),]
data_b <- b %>% pivot_longer(c(3:5),names_to = "model", values_to = "values")
cor_total <- data_b
cor_total$tool <- "TOTAL"
boxplot(values~model,data=cor_total)

#sedi
b <- left_join(ds[,c(1,10,11)], gl[,c(1,10,11)], c("lead", "sedi_months"))
colnames(b)[3:4] <- c("Downscaled", "Global Full")
b <- left_join(b, gl_r[,c(1,10,11)], c("lead", "sedi_months"))
colnames(b)[5] <- c("Global Reduced")
colnames(b)[2] <- c("month")
b <- b[!duplicated(b),]
data_b <- b %>% pivot_longer(c(3:5),names_to = "model", values_to = "values")
sedi_total <- data_b
sedi_total$tool <- "TOTAL"
boxplot(values~model,data=sedi_total)

#ggplot solution
acc$metric <- "Forecast Accuracy"
cor$metric <- "Correlation"
sedi$metric <- "SEDI"
hci <- rbind(cor, acc, sedi)

acc_total$metric <- "Forecast Accuracy"
cor_total$metric <- "Correlation"
sedi_total$metric <- "SEDI"
total <- rbind(cor_total, acc_total, sedi_total)

total <- total[,c(2,1,3,4,5,6)]

data_all <- rbind(hci, total)

g1 <- ggplot(data_all,aes(x=model, y=values))+
  geom_boxplot()+
  geom_jitter(aes(x=model, y=values), col="grey", size=1, alpha=0.9, width=0.1) +
  facet_wrap(~tool, nrow=2, scales="free")+
  ylab("Skill Assessment Value")+xlab("")+
  theme_bw()+
  facet_grid(tool~metric, scales="free")

tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/Fig4_skill_boxplots.tiff', units="in",width=9, height=5, res=330)
plot(g1)
dev.off()



#POSTER PLOT
# 
# g1_poster <- ggplot(hci[hci$metric=="Correlation",],aes(x=model, y=values))+
#   geom_boxplot()+
#   geom_jitter(aes(x=model, y=values), col="grey", size=1, alpha=0.9, width=0.1) +
#   facet_wrap(~tool, nrow=2, scales="free")+
#   ylab("Correlation")+xlab("")+
#   theme_bw()+
#   theme(rect = element_rect(fill = "transparent"))
# 
# tiff('~/Dropbox/Conference and Workshops/AMSA 2023/dscale_plot.tiff', units="in",width=4, height=3, res=330, bg="transparent")
# plot(g1_poster)
# dev.off()

#-----Figure S1: Comparative & reduced scenario of HCI--------
##skill maps
metric_summary <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_3b/global_forecasts_metric_summary.rds')
metric_summary_init1 <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_5/downscaled_forecast_metric_summary_init1.rds')
metric_summary_init7 <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_5/downscaled_forecast_metric_summary_init7.rds')

metric_summary_init1 <- na.omit(metric_summary_init1)
metric_summary_init7 <- na.omit(metric_summary_init7)
ds <- rbind(metric_summary_init1, metric_summary_init7)


# #Correlation plots
# cor <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=cor))+
#   geom_tile(color="white")+
#   scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
#                        midpoint = 0.36, limit = c(0,1), space = "Lab", 
#                        name="") +
#   scale_x_continuous(breaks=1:12,expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(x="Month", y="Lead time (months)")+
#   ggtitle("Correlation: global model")+
#   theme(legend.position = "right",
#         legend.key.height= unit(0.5, 'cm'),
#         legend.key.width= unit(0.2, 'cm'),
#         plot.title = element_text(size = 10))
#ACC
acc <- ggplot(data = metric_summary, aes(x=month,y=lead,fill=acc))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(0,0.5,1)), limit = c(0,1), space = "Lab",
                       name="") +
  scale_x_continuous(breaks = seq_along(month.abb), labels = month.abb, expand = c(0, 0))+
  scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
  labs(x="", y="")+
  ggtitle("Global Full Ensemble")+
  theme(legend.position = "right",
        legend.key.height= unit(0.6, 'cm'),
        legend.key.width= unit(0.3, 'cm'),
        plot.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45))
# (0.5*0.5) + (0.5*0.5) = 0.5 #based on median (50% quantile) being roughly 0.5

#DOWNSCALED Correlation plots

#ACC
ds_acc <- ggplot(data = ds, aes(x=month,y=lead,fill=acc))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(0,0.5,1)), limit = c(0,1), space = "Lab",
                       name="") +
  scale_x_continuous(breaks = seq_along(month.abb), labels = month.abb, expand = c(0, 0))+
  scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
  # labs(x="Month", y="Lead time (months)")+
  labs(x="", y="")+
  ggtitle("Downscaled Ensemble")+
  theme(legend.position = "right",
        legend.key.height= unit(0.6, 'cm'),
        legend.key.width= unit(0.3, 'cm'),
        plot.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45))

#Difference plots:
ds$model <- "downscaled"
metric_summary$model <- "global"
all <- rbind(ds, metric_summary)
# cor_long <- all %>% pivot_wider(id_cols=1:2,names_from = model,values_from=cor)
acc_long <- all %>% pivot_wider(id_cols=1:2,names_from = model,values_from=acc)
# cor_long$diff <- cor_long$global-cor_long$downscaled
acc_long$diff <- acc_long$global-acc_long$downscaled
# cor_long <- na.omit(cor_long)
acc_long <- na.omit(acc_long)

g4 <- ggplot(acc_long, aes(x=month, y=lead, fill=diff))+
  geom_tile(color="white")+
  geom_text(aes(label = round(diff, 2)), size=2) +
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(-0.33,0,0.33)), limit = c(-0.33,0.33), space = "Lab",
                       name="") +
  scale_x_continuous(breaks = seq_along(month.abb), labels = month.abb, expand = c(0, 0))+
  scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
  ggtitle("Difference: global full - downscaled")+
  labs(x="", y="")+
  theme(legend.position = "right",
        legend.key.height= unit(0.6, 'cm'),
        legend.key.width= unit(0.3, 'cm'),
        plot.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45))

#NOw get REDUCED scenario
##skill maps
metric_summary_r <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_3c/global_forecasts_metric_summary_only3.rds')
metric_summary_init1 <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_5/downscaled_forecast_metric_summary_init1.rds')
metric_summary_init7 <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/HCI/step_5/downscaled_forecast_metric_summary_init7.rds')
metric_summary_init1 <- na.omit(metric_summary_init1)
metric_summary_init7 <- na.omit(metric_summary_init7)
ds <- rbind(metric_summary_init1, metric_summary_init7)

#ACC
acc_r <- ggplot(data = metric_summary_r, aes(x=month,y=lead,fill=acc))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(0,0.5,1)), limit = c(0,1), space = "Lab",
                       name="") +
  scale_x_continuous(breaks = seq_along(month.abb), labels = month.abb, expand = c(0, 0))+
  scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
  labs(x="", y="")+
    ggtitle("Global Reduced Ensemble")+
  theme(legend.position = "right",
        legend.key.height= unit(0.6, 'cm'),
        legend.key.width= unit(0.3, 'cm'),
        plot.title = element_text(size = 10),
        axis.text.x = element_text(angle = 45))
# (0.5*0.5) + (0.5*0.5) = 0.5 #based on median (50% quantile) being roughly 0.5

#Difference plots:
ds$model <- "downscaled"
metric_summary_r$model <- "global"
all_r <- rbind(ds, metric_summary_r)
# cor_long <- all %>% pivot_wider(id_cols=1:2,names_from = model,values_from=cor)
acc_long_r <- all_r %>% pivot_wider(id_cols=1:2,names_from = model,values_from=acc)
# cor_long$diff <- cor_long$global-cor_long$downscaled
acc_long_r$diff <- acc_long_r$global-acc_long_r$downscaled
# cor_long <- na.omit(cor_long)
acc_long_r <- na.omit(acc_long_r)

g4_r <- ggplot(acc_long_r, aes(x=month, y=lead, fill=diff))+
  geom_tile(color="white")+
  geom_text(aes(label = round(diff, 2)), size=2) +
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(-0.33,0,0.33)), limit = c(-0.33,0.33), space = "Lab",
                       name="") +
  scale_x_continuous(breaks = seq_along(month.abb), labels = month.abb, expand = c(0, 0))+
  scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
  ggtitle("Difference: global reduced - downscaled")+
  labs(x="", y="")+
  theme(legend.position = "right",
        legend.key.height= unit(0.6, 'cm'),
        legend.key.width= unit(0.3, 'cm'),
        plot.title = element_text(size = 9),
        axis.text.x = element_text(angle = 45))


hlay <- rbind(c(1,2),
              c(3,NA),
              c(4,5))
tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/FigS2_HCI_skill_comparitive_reduced.tiff', units="in", res=300, width=7, height=8)
grid.arrange(grobs=list(acc, acc_r,
                        ds_acc,
                        g4,g4_r),
             left="Lead time (months)", bottom="Forecast Month", 
             layout_matrix=hlay)
dev.off()


#-----Figure S2: comparative & reduced scenario of TOTAL - Forecast Accuracy------
##skill maps for comparative scenario (i.e. 73 ensemble members)
gl <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_8b/metric_summary_global_1980-2010_applesscenario.rds')
ds <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_9b/metric_downscaled.rds')
ds <- na.omit(ds)

#ACC
gl$acc_months <- as.factor(gl$acc_months)
gl$acc_months <- factor(gl$acc_months, levels=rev(levels(gl$acc_months)))
gl$lead <- gl$lead+6
acc <-   ggplot(data = gl, aes(x=acc_months,y=lead,fill=acc_values))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(0,0.62,1)), limit = c(0,1), space = "Lab",
                       name="") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
  labs(x="", y="")+
  ggtitle("Global Full Ensemble")+
  theme(legend.position = "right",
        legend.key.height= unit(0.6, 'cm'),
        legend.key.width= unit(0.3, 'cm'),
        plot.title = element_text(size = 10))
# (0.5*0.5) + (0.5*0.5) = 0.5 #based on median (50% quantile) being roughly 0.5


#DOWNSCALED ACC plots
ds$acc_months <- as.factor(ds$acc_months)
ds$acc_months <- factor(ds$acc_months, levels=rev(levels(ds$acc_months)))
ds$lead <- ds$lead+6
ds_acc <- ggplot(data = ds, aes(x=acc_months,y=lead,fill=acc_values))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(0,0.71,1)), limit = c(0,1), space = "Lab",
                       name="")+
  scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
  # labs(x="Closure Month", y="Lead time (months)")+
  labs(x="", y="")+
  ggtitle("Downscaled Ensemble")+
    theme(legend.position = "right",
          legend.key.height= unit(0.6, 'cm'),
          legend.key.width= unit(0.3, 'cm'),
          plot.title = element_text(size = 10))
# (0.5*0.5) + (0.5*0.5) = 0.5 #based on median (50% quantile) being roughly 0.5

#Difference plots:
ds$model <- "downscaled"
gl$model <- "global"
all <- rbind(ds, gl)

acc_long <- all[,c(1,8,9,14)]  %>% distinct() %>% pivot_wider(names_from = model,values_from=acc_values)
acc_long$diff <- acc_long$global-acc_long$downscaled
acc_long <- na.omit(acc_long)

acc_long$acc_months <- as.factor(acc_long$acc_months)
acc_long$acc_months <- factor(acc_long$acc_months, levels=rev(levels(acc_long$acc_months)))
g4 <- ggplot(acc_long, aes(x=acc_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(-0.14,0,0.14)), limit = c(-0.14,0.14), space = "Lab",
                       name="")+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
  ggtitle("Difference: global full - downscaled")+
  labs(x="", y="")+
  theme(legend.position = "right",
        legend.key.height= unit(0.6, 'cm'),
        legend.key.width= unit(0.3, 'cm'),
        plot.title = element_text(size = 10))

#REDUCED SCENARIO
gl_r <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_8c/metric_summary_global_1980-2010_only3_scenario.rds')
ds_r <- readRDS('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/TOTAL/step_9b/metric_downscaled.rds')
ds_r <- na.omit(ds_r)

#ACC
gl_r$acc_months <- as.factor(gl_r$acc_months)
gl_r$acc_months <- factor(gl$acc_months, levels=rev(levels(gl_r$acc_months)))
gl_r$lead <- gl_r$lead+6
acc_r <- ggplot(data = gl_r, aes(x=acc_months,y=lead,fill=acc_values))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(0,0.71,1)), limit = c(0,1), space = "Lab",
                       name="")+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
  labs(x="", y="")+
  ggtitle("Global Reduced Ensemble")+
  theme(legend.position = "right",
        legend.key.height= unit(0.6, 'cm'),
        legend.key.width= unit(0.3, 'cm'),
        plot.title = element_text(size = 10))
# (0.5*0.5) + (0.5*0.5) = 0.5 #based on median (50% quantile) being roughly 0.5

#Difference plots:
gl_r$model <- "global"
all_r <- rbind(ds, gl_r)

acc_long_r <- all_r[,c(1,8,9,14)]  %>% distinct() %>% pivot_wider(names_from = model,values_from=acc_values)
acc_long_r$diff <- acc_long_r$global-acc_long_r$downscaled
acc_long_r <- na.omit(acc_long_r)

acc_long_r$acc_months <- as.factor(acc_long_r$acc_months)
acc_long_r$acc_months <- factor(acc_long_r$acc_months, levels=rev(levels(acc_long_r$acc_months)))
g4_r <- ggplot(acc_long_r, aes(x=acc_months, y=lead, fill=diff))+
  geom_tile()+
  geom_text(aes(label = round(diff, 2))) +
  scale_fill_gradientn(colors = cmocean("balance")(256),
                       values=rescale(c(-0.14,0,0.14)), limit = c(-0.14,0.14), space = "Lab",
                       name="")+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(breaks=0.5:11.5,expand = c(0, 0)) +
  ggtitle("Difference: global reduced - downscaled")+
  labs(x="", y="")+
  # labs(x="Closure Month", y="Lead time (months)")+
  theme(legend.position = "right",
        legend.key.height= unit(0.6, 'cm'),
        legend.key.width= unit(0.3, 'cm'),
        plot.title = element_text(size = 10))



hlay <- rbind(c(1,2),
              c(3,NA),
              c(4,5))

tiff('~/Dropbox/PROJECTS/MAPP/MAPP_EcoForecast/manuscript/figures/FigS3_TOTAL_skill_reduced.tiff', units="in", res=300, width=7, height=8)
grid.arrange(grobs=list(acc, acc_r,
          ds_acc,
          g4,g4_r),
          left="Lead time (months)", bottom="Closure Month", 
          layout_matrix=hlay)
dev.off()










