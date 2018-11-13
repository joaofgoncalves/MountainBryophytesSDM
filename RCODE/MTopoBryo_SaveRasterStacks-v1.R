

#projNames

CoordRefSyst <- crs(current)

for(i in 1:length(projNames)){
  
  tmp <- projectRaster(from = get(projNames[i]), res = 0.008333334, crs=CoordRefSyst, method = 'ngb')
  
  writeRaster(tmp, filename = paste("./DATA/RASTER/_VARS/",projNames[i],".tif",sep=""), overwrite=TRUE)
  
  print(projNames[i])
  
}
