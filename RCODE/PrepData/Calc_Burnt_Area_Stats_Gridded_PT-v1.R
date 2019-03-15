
library(sf)
library(fasterize)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(dplyr)

aa <- read_sf("D:/MyDocs/Projects/ERPAM_FloraData/DATA/VECTOR/Areas_ardidas_ICNF_1990_2017/AreaArdidaICNF_1990_2017.shp")

PTbounds <- readOGR("D:/MyDocs/Projects/EcosystemStability_RS_PT/DATA/_TMP","PTbounds_PAs_v1")

# r <- raster(ext=extent(PTbounds), resolution=50, crs=crs(PTbounds))
# 
# i <- 0
# yrs <- 2000:2017
# 
# for(yr in yrs){
#   
#   i<-i+1
#   
#   cat("[",i,"/",length(yrs),"] Rasterizing year:",yr,".....")
#   
#   # Subset vector layer to each year and merge it
#   #tmp <- gUnaryUnion(shp[shp$YEAR == yr,])
#   tmp <- aa[aa$YEAR == yr, ]
#   
#   # Rasterize object | 0 = not burnt / 1 = burned
#   rst <- fasterize(tmp, r, field = NULL, background=0)
#   
#   # if(i==1)
#   #   rstBurntFreq <- rst
#   # else
#   #   rstBurntFreq <- stack(rstBurntFreq, rst)
#   
#   writeRaster(rst, filename = paste("D:/MyDocs/temp/tmp_HH/ba/burntAreaPT_yr",yr,".tif",sep=""))
#   
#   
#   cat("done.\n\n")
#   
#   
# }
# 
# fl <- list.files("D:/MyDocs/temp/tmp_HH/ba",pattern=".tif$",full.names = TRUE)
# rstBurntFreq <- stack(fl)
# 
# names(rstBurntFreq) <- paste("yr_",yrs,sep="")
# 
# freqSum_00_17 <- calc(rstBurntFreq, fun = sum)
# writeRaster(freqSum_00_17, filename = "D:/MyDocs/temp/tmp_HH/recAreaArdida_2000_2017_PT.tif")



grid_pt <- read_sf("D:/MyDocs/temp/tmp_HH/grid1km/grid1km_PT.shp", stringsAsFactors = FALSE)

i <- 0
yrs <- 2000:2017

baGrid <- data.frame(Id = grid_pt$Id)


for(yr in yrs){
  
  i<-i+1
  
  cat("[",i,"/",length(yrs),"] Calculating year:",yr,".....")
  
  # Subset vector layer to each year and merge it
  tmp <- st_buffer(aa[aa$YEAR == yr, ],dist = 0)
  
  ints <- st_intersects(grid_pt, tmp, sparse = FALSE, prepared = TRUE)
  
  int <- st_intersection(tmp, grid_pt[apply(ints,1,sum)>0,])

  bas <- cbind(Id = int$Id, area = st_area(int)) %>% 
    as.data.frame %>% 
    group_by(Id) %>% 
    summarise(tarea=sum(area))
  
  baGrid <- merge(baGrid,bas,by="Id",all=TRUE)
  
  cat("done.\n\n")
  
}

baGrid[is.na(baGrid)] <- 0
colnames(baGrid) <- c("Id", paste("ba_yr",yrs,sep=""))

ba_avg <- data.frame(Id = baGrid[,1], ba_mean = (apply(baGrid[,-1],1,mean)/1E6)*100)
write.csv(ba_avg, "D:/MyDocs/temp/tmp_HH/ba_avg.csv",row.names = FALSE)

