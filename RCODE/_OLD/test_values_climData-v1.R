

library(raster)
library(magrittr)


dirList <- list.dirs("./DATA/RASTER/Worldclim_1km")[-1]
projList <- c()


fl.topoVars <- list.files("./DATA/RASTER/TopoVars", full.names = TRUE, pattern = ".tif$")[c(1,2,3,5)]


i <- 0
for(fld in dirList){
  
  cat("-> Loading data for:",fld,"..... \n")  
  i<-i+1
  projList[i] <- basename(fld)
  fl.climVars <- list.files(fld, pattern = ".tif$", full.names = TRUE)
  
  print(fl.climVars)
  
  #try(assign(x = projList[i], value = stack(c(fl.topoVars, fl.climVars))))
  try(assign(x = projList[i], value = stack(c(fl.climVars))))
  
  cat("done.\n\n")
}








pt <- SpatialPoints(data.frame(x=c(0.912479, 
                                   -1.778232,
                                   -5.473310),
                               y=c(42.561549,
                                   40.453618,
                                   40.242825)))

# Current
extract(a2000, pt) %>% round(2)

# he45
extract(he45bi50, pt) %>% round(2)
extract(he45bi70, pt) %>% round(2)

# he85
extract(he85bi50, pt) %>% round(2)
extract(he85bi70, pt) %>% round(2)

# mp45
extract(mp45bi50, pt) %>% round(2)
extract(mp45bi70, pt) %>% round(2)

#mp85
extract(mp85bi50, pt) %>% round(2)
extract(mp85bi70, pt) %>% round(2)






compareRaster(a2000[[5]], he45bi50[[5]], values = TRUE)
compareRaster(a2000[[6]], he45bi50[[6]], values = TRUE)
compareRaster(a2000[[7]], he45bi50[[7]], values = TRUE)
compareRaster(a2000[[8]], he45bi50[[8]], values = TRUE)
compareRaster(a2000[[9]], he45bi50[[9]], values = TRUE)


compareRaster(a2000[[5]], he45bi70[[5]], values = TRUE)
compareRaster(a2000[[6]], he45bi70[[6]], values = TRUE)
compareRaster(a2000[[7]], he45bi70[[7]], values = TRUE)
compareRaster(a2000[[8]], he45bi70[[8]], values = TRUE)
compareRaster(a2000[[9]], he45bi70[[9]], values = TRUE)


compareRaster(a2000[[5]], he85bi50[[5]], values = TRUE)
compareRaster(a2000[[6]], he85bi50[[6]], values = TRUE)
compareRaster(a2000[[7]], he85bi50[[7]], values = TRUE)
compareRaster(a2000[[8]], he85bi50[[8]], values = TRUE)
compareRaster(a2000[[9]], he85bi50[[9]], values = TRUE)










