
library(readxl)
library(raster)
library(dplyr)
library(sf)
library(lwgeom)


setwd("D:/MyDocs/Dropbox/Modelling_mountaintop_bryophytes_climate_change/MountainBryophytesSDM/OUT/MODS/BIN_PREDS")

fl <- list.files(pattern=".tif$")

spData <- read_sf("D:/MyDocs/Dropbox/Modelling_mountaintop_bryophytes_climate_change/MountainBryophytesSDM/DATA/VECTOR/Bryophyte_dataset/And_Gri_Rac_PI_all_2.shp")

spNames <- sort(unique(spData$Cod_esp))

rstBinMods <- raster::stack(fl[1])


N <- nlayers(rstBinMods)

wcentr <- matrix(NA,N,2)
centr <- matrix(NA,N,2)

for(i in 1:N){
  
  rstBin <- rstBinMods[[i]]
  
  rstBin[rstBin == 0] <- NA
  
  #plot(rstBin)
  
  rstPol <- rasterToPolygons(rstBin, n = 8)
  rstPol_sf <- st_as_sf(rstPol) 
  rstPol_sf <- rstPol_sf %>% 
    st_transform(crs="+proj=utm +zone=30 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") %>% 
    mutate(area_m2 = st_area(rstPol_sf))
  
  sfpts <- st_cast(rstPol_sf, 'POINT')
  sfCoords <- st_coordinates(sfpts)
  
  wcentr[i,] <- apply(cbind(sfCoords[,"X"]*sfpts$area_m2,sfCoords[,"Y"]*sfpts$area_m2), 2, sum, na.rm=TRUE)/sum(sfpts$area_m2)
  centr[i,] <- apply(sfCoords, 2, mean)
  
  cat("Finished",i,"/",N,"..\nSp. name:",spNames[i],"\n\n")

}

wcentr_ <- data.frame(spName = spNames, wcentr, stringsAsFactors = FALSE)
centr_ <- data.frame(spName = spNames, centr, stringsAsFactors = FALSE)

write.csv(wcentr_, file = "../../WeightedSpatialCenters-v1.csv", row.names = FALSE)
write.csv(centr_, file = "../../SpatialCenters-v1.csv", row.names = FALSE)

d1 <- sqrt(apply((apply(wcentr,2,mean) - wcentr)^2, 1, sum))
d2 <- sqrt(apply((apply(centr,2,mean) - centr)^2, 1, sum))


hist(log(d1))
hist(log(d2))


spatialPositionDF <- data.frame(spName = spNames, 
                                spPos=d2, spPosLog = log(d2),
                                wspPos=d1, wspPosLog = log(d1),
                                stringsAsFactors = FALSE)

write.csv(spatialPositionDF, file = "../../SpatialPositionIndices-v1.csv", row.names = FALSE)




