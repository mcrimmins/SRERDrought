##### Download PRISM Climate and calculate drought indices
# MAC 08/20/24

# load libraries
library(sf)
library(RCurl)
library(jsonlite)
library(raster)
library(SCI)
library(terra)
library(SPEI)
library(tidyr)
library(dplyr)

# import shapefile
shp<- st_read(dsn = paste0(getwd(),"/shapefiles/SLUD"))
plot(st_geometry(shp)) # plot to check

# get bbox
bbox<-as.numeric(st_bbox(shp))
# convert to ACISbbox
ACISbbox<-paste(as.character(bbox), collapse=",")
#ACISbbox<- "-111.45093,31.33234,-108.92062,33.00649"

# set download parameters -- period of record for PRISM 
startYr<-1895
endYr<-2023
endMo<-12
endDay<-1

# create current date and generate dates -- keep with PRISM date
dateRangeStart = paste0(startYr,"-01-01")
dateRangeEnd = paste0(endYr,"-",endMo,"-",endDay)
allDates<-seq(as.Date(dateRangeStart), as.Date(dateRangeEnd),by="month")

# RCC-ACIS query
jsonQuery=paste0('{"bbox":"',ACISbbox,'","sdate":"',dateRangeStart,'","edate":"',dateRangeEnd,'","grid":"21","elems":[{"name":"mly_pcpn","units":"mm"},{"name":"mly_avgt","units":"degreeC"},{"name":"mly_mint","units":"degreeC"},{"name":"mly_maxt","units":"degreeC"}],"meta":"ll,elev","output":"json"}') # or uid
# submit
out<-postForm("http://data.rcc-acis.org/GridData",
              .opts = list(postfields = jsonQuery,
                           httpheader = c('Content-Type' = 'application/json', Accept = 'application/json')))
# convert from JSON
out<-fromJSON(out)

#### loop through vars, process into gridstacks and place in list ####
gridList<-list()
for(k in 2:5){
  # convert to list of matrices, flipud with PRISM
  matrixList <- vector("list",length(out$data))
  for(i in 1:length(out$data)){
    matrixList[[i]]<-apply(t(out$data[[i]][[k]]),1,rev)
  }
  # read into raster stack
  rasterList<-lapply(matrixList, raster)
  gridStack<-stack(rasterList)
  gridExtent<-extent(min(out$meta$lon), max(out$meta$lon), min(out$meta$lat), max(out$meta$lat))
  gridStack<-setExtent(gridStack, gridExtent, keepres=FALSE, snap=FALSE)
  names(gridStack)<-allDates
  # put into list
  gridList[[k-1]]<-gridStack
}
#####

##### layers into time series
# process precip
pptGrid<-gridList[[1]]
# set 0 and neg to NA
pptGrid[pptGrid < 0] <- NA
# convert to terra format
pptGrid<-rast(pptGrid)
# get mean by layer
pptTS<-global(pptGrid, "mean", na.rm=TRUE)

# process tmean
tmeanGrid<-gridList[[2]]
# set 0 and neg to NA
tmeanGrid[tmeanGrid <= -999] <- NA
# convert to terra format
tmeanGrid<-rast(tmeanGrid)
# get mean by layer
tmeanTS<-global(tmeanGrid, "mean", na.rm=TRUE)

# process tmin
tminGrid<-gridList[[3]]
# set 0 and neg to NA
tminGrid[tminGrid <= -999] <- NA
# convert to terra format
tminGrid<-rast(tminGrid)
# get mean by layer
tminTS<-global(tminGrid, "mean", na.rm=TRUE)

# process tmax
tmaxGrid<-gridList[[4]]
# set 0 and neg to NA
tmaxGrid[tmaxGrid <= -999] <- NA
# convert to terra format
tmaxGrid<-rast(tmaxGrid)
# get mean by layer
tmaxTS<-global(tmaxGrid, "mean", na.rm=TRUE)

# combine into dataframe
climateData<-cbind.data.frame(allDates,as.numeric(format(allDates,"%m")),tmeanTS$mean,tmaxTS$mean,tminTS$mean,pptTS$mean)
colnames(climateData)<-c("date","month","tmean","tmax","tmin","precip")
#####

# calculate SPI using SCI package
spiList<-list()
# loop through SPI timescales
for(i in 1:120){
  spiList[[i]]<-transformSCI(climateData$precip,first.mon=1,
                                  obj=fitSCI(climateData$precip,first.mon=1,distr="gamma",
                                                  time.scale=i,p0=TRUE),sci.limit=3)
}
#  convert to dataframe
spiDF<-do.call(cbind,spiList)
  spiDF<-cbind.data.frame(climateData$date,spiDF)
  colnames(spiDF)[1]<-c("date")
  colnames(spiDF)[2:ncol(spiDF)]<-paste0("SPI-",colnames(spiDF)[2:ncol(spiDF)])
  
# calculate SPEI using SPEI::Hargreaves and SCI gev distribution
climateData$PET<-hargreaves(climateData$tmin,climateData$tmax,lat=(bbox[2]+bbox[4])/2)
climateData$waterBalance<-climateData$precip-climateData$PET
speiList<-list()
# loop through SPEI timescales
# library(evd) ## add in for gev distribution
for(i in 1:120){
 speiList[[i]]<-transformSCI(climateData$waterBalance,first.mon=1,
                            obj=fitSCI(climateData$waterBalance,first.mon=1,distr="genlog",
                                       time.scale=i,p0=FALSE,start.fun.fix = FALSE),sci.limit=3)
  #speiList[[i]]<-spei(climateData$waterBalance, i)$fitted
}
#  convert to dataframe
speiDF<-do.call(cbind,speiList)
speiDF<-cbind.data.frame(climateData$date,speiDF)
colnames(speiDF)[1]<-c("date")
colnames(speiDF)[2:ncol(speiDF)]<-paste0("SPEI-",colnames(speiDF)[2:ncol(speiDF)])

#####
# normality tests of all SPI and SPEI calculations
# following Wu et al. 2007 https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1002/joc.1371

# spi tests
spiDF$month<-as.numeric(format(spiDF$date,"%m"))
spiDFlong<-spiDF %>% pivot_longer(cols=`SPI-1`:`SPI-120`,
                                  names_to = 'timescale',
                                  values_to = 'values')
spiNormStats<-spiDFlong %>% group_by(month,timescale) %>%
                            summarize(sh_pval=shapiro.test(values)$p.value,
                                      sh_stat=shapiro.test(values)$statistic,
                                      medianVal=median(values, na.rm=TRUE))
spiNormStats$pval_flag<-ifelse(spiNormStats$sh_pval<0.1, 1,0)
spiNormStats$stat_flag<-ifelse(spiNormStats$sh_stat<0.96, 1,0)
spiNormStats$med_flag<-ifelse(abs(spiNormStats$medianVal)>0.05, 1,0)
spiNormStats$flagSum<-spiNormStats$pval_flag+spiNormStats$stat_flag+spiNormStats$med_flag
# look at flagged distributions
table(spiNormStats$flagSum)
# hist(subset(spiDF[,c("SPI-2","month")], month==2)$`SPI-2`)
# hist(subset(spiDF[,c("SPI-1","month")], month==4)$`SPI-1`)
#   hist(subset(climateData, month==4)$precip)

# spei tests
speiDF$month<-as.numeric(format(speiDF$date,"%m"))
speiDFlong<-speiDF %>% pivot_longer(cols=`SPEI-1`:`SPEI-120`,
                                  names_to = 'timescale',
                                  values_to = 'values')
speiNormStats<-speiDFlong %>% 
  drop_na(values) %>%
  group_by(month,timescale) %>%
                            summarize(sh_pval=shapiro.test(values)$p.value,
                                      sh_stat=shapiro.test(values)$statistic,
                                      medianVal=median(values, na.rm=TRUE))
speiNormStats$pval_flag<-ifelse(speiNormStats$sh_pval<0.1, 1,0)
speiNormStats$stat_flag<-ifelse(speiNormStats$sh_stat<0.96, 1,0)
speiNormStats$med_flag<-ifelse(abs(speiNormStats$medianVal)>0.05, 1,0)
speiNormStats$flagSum<-speiNormStats$pval_flag+speiNormStats$stat_flag+speiNormStats$med_flag
# look at flagged observations
table(speiNormStats$flagSum)
# hist(subset(speiDF[,c("SPEI-120","month")], month==2)$`SPEI-120`)
# hist(subset(spiDF[,c("SPI-1","month")], month==4)$`SPI-1`)
# hist(subset(climateData, month==4)$precip)

# write out dataframes to Excel spreadsheet
# library(xlsx)
# write.xlsx(spiDF, file = "SRER_DroughtIndices.xlsx",
#            sheetName = "SPI", append = FALSE)
# # Add SPEI
# write.xlsx(speiDF, file = "SRER_DroughtIndices.xlsx", 
#            sheetName="SPEI", append=TRUE)
# # Add SPI diags
# write.xlsx(spiNormStats, file = "SRER_DroughtIndices_Diags.xlsx",
#            sheetName="SPI_diagnostics", append=FALSE)
# # Add SPEI diags
# write.xlsx(speiNormStats, file = "SRER_DroughtIndices_Diags.xlsx",
#            sheetName="SPEI_diagnostics", append=TRUE)

# write out to csv files
write.csv(climateData, file = "SRER_Climate.csv")
write.csv(spiDF, file = "SPI.csv")
write.csv(speiDF, file = "SPEI.csv")
write.csv(spiNormStats, file = "SPI_diag.csv")
write.csv(speiNormStats, file = "SPEI_diag.csv")


  