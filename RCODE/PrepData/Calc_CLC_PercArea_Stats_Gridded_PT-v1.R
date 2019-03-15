
library(sf)
library(fasterize)
library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(dplyr)

aa <- read_sf("D:/MyDocs/temp/colabHH/clc2012/clc2012_WGS84UTM29N.shp")

#PTbounds <- readOGR("D:/MyDocs/Projects/EcosystemStability_RS_PT/DATA/_TMP","PTbounds_PAs_v1")

grid_pt <- read_sf("D:/MyDocs/Dropbox/Modelling_mountaintop_bryophytes_climate_change/MountainBryophytesSDM/DATA/VECTOR/grid1km/grid1km_PT.shp", stringsAsFactors = FALSE)

i <- 0
clcClasses <- c("312","313","331","333","321","322","323","324")

baGrid <- data.frame(Id = grid_pt$Id)


for(clcClass in clcClasses){
  
  i<-i+1
  
  cat("[",i,"/",length(clcClasses),"] Calculating class:",clcClass,".....")
  
  # Subset vector layer to each year and merge it
  tmp <- st_buffer(aa[aa$code_12 == clcClass, ],dist = 0)
  
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
colnames(baGrid) <- c("Id", paste("clc_",clcClasses,sep=""))

scrHerbVeg <- apply(baGrid[,6:9], 1, sum)

clcGrid <- data.frame(baGrid, clc_shv = scrHerbVeg)

clcGrid <- data.frame(Id=clcGrid$Id, (clcGrid[,-1]/1E6)*100)


write.csv(clcGrid, "D:/MyDocs/temp/colabHH/clcGrid.csv",row.names = FALSE)



