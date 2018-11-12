
library(dplyr)
library(raster)
library(biomod2)
library(biomod2plus)
library(readxl)
library(sp)
library(rgdal)
library(magrittr)
library(tools)
library(stringr)



#setwd("D:/MyDocs/Dropbox/Modelling_mountaintop_bryophytes_climate_change/MountainBryophytesSDM")
setwd("/home/cibio/myfiles/MountainBryophytesSDM")

## -------------------------------------------------------------------------------------- ##
## LOAD INPUT PRESENCE POINTS ----
## -------------------------------------------------------------------------------------- ##

# Read data in shapefile format
spBryoData <- readOGR("./DATA/VECTOR/Bryophyte_dataset","And_Gri_Rac_PI_all_2")

# Make a data.frame object with xy coordinates
spBryoDF <- data.frame(coordinates(spBryoData),spBryoData@data)
colnames(spBryoDF)[1:2] <- c("lon","lat")

# Species codes to model
allSpNames <- unique(spBryoData$Cod_esp)

## -------------------------------------------------------------------------------------- ##
## LOAD RASTER DATA  ----
## -------------------------------------------------------------------------------------- ##

# Data paths to raster variables
#climVarsPaths <- list.files("./DATA/RASTER/WorldClim_1km/a2000", pattern=".tif$", full.names = TRUE)
# climVarsPaths <- list.files("/home/cibio/myfiles/MountainBryophytesSDM/DATA/RASTER/Clim/a2000", pattern=".tif$", full.names = TRUE)
# topoVarsPaths <- list.files("/home/cibio/myfiles/MountainBryophytesSDM/DATA/RASTER/TopoVars", pattern=".tif$", full.names = TRUE)[c(2,5)]
# soilVarsPaths <- list.files("/home/cibio/myfiles/MountainBryophytesSDM/DATA/RASTER/Soil", pattern=".tif$", full.names = TRUE)[c(8)]
# 
# # Load raster stack with climatic and topographic data
# current <- stack(c(climVarsPaths,topoVarsPaths,soilVarsPaths))
# names(current) <- c(paste("BIO",c(11,17,19),sep=""),"ASPBR","TOPRI","SOIPH")
# 
# # Directories 
# climProjDirs <- list.dirs("/home/cibio/myfiles/MountainBryophytesSDM/DATA/RASTER/Clim", recursive = FALSE)[-c(1,seq(3,9,by=2))]
# 
# climProjNames <- basename(climProjDirs)
# 
# 
# he45bi50 <- stack(stack("/home/cibio/myfiles/MountainBryophytesSDM/DATA/RASTER/Clim/he45bi50_proj_BIO_11_17_19.tif"),
#                   stack(topoVarsPaths,soilVarsPaths))
# names(he45bi50) <- c(paste("BIO",c(11,17,19),sep=""),"ASPBR","TOPRI","SOIPH")
# 
# 
# he85bi50 <- stack(stack("/home/cibio/myfiles/MountainBryophytesSDM/DATA/RASTER/Clim/he85bi50_proj_BIO_11_17_19.tif"),
#                   stack(topoVarsPaths,soilVarsPaths))
# names(he85bi50) <- c(paste("BIO",c(11,17,19),sep=""),"ASPBR","TOPRI","SOIPH")
# 
# 
# mp45bi50 <- stack(stack("/home/cibio/myfiles/MountainBryophytesSDM/DATA/RASTER/Clim/mp45bi50_proj_BIO_11_17_19.tif"),
#                   stack(topoVarsPaths,soilVarsPaths))
# names(mp45bi50) <- c(paste("BIO",c(11,17,19),sep=""),"ASPBR","TOPRI","SOIPH")
# 
# 
# mp85bi50 <- stack(stack("/home/cibio/myfiles/MountainBryophytesSDM/DATA/RASTER/Clim/mp85bi50_proj_BIO_11_17_19.tif"),
#                   stack(topoVarsPaths,soilVarsPaths))
# names(mp85bi50) <- c(paste("BIO",c(11,17,19),sep=""),"ASPBR","TOPRI","SOIPH")



fl <- list.files("./DATA/RASTER/_VARS", pattern=".tif$", full.names = TRUE)
projNames <- gsub(".tif","",basename(fl))

i <- 0
for(fn in fl){
  i<-i+1
  assign(projNames[i], stack(fn))
  assign(projNames[i],`names<-`(get(projNames[i]), 
                                c(paste("BIO",c(11,17,19),sep=""),"ASPBR","TOPRI","SOIPH")))
  print(projNames[i])
}



## -------------------------------------------------------------------------------------- ##
## PARAMETERS ----
## -------------------------------------------------------------------------------------- ##

projNames <- c("current", climProjNames)




sp <- "GRIMER"

setwd("/mnt/dados/datasets/jg/OUT/MODS")

load("/mnt/dados/datasets/jg/OUT/MODS/GRIMER/GRIMER.1541781256.models.out")

myBiomodModelOut <- GRIMER.1541781256.models.out




# Ensemble all partial models
selMods <- twoStepBestModelSelection(myBiomodModelOut, 
                                     evalMetric = "TSS", 
                                     nrBestAlgos = 7, 
                                     bestAlgoFun = stats::median, 
                                     topFraction = 0.05)

myBiomodEM <- try(BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                          chosen.models = selMods,
                                          em.by = 'all',
                                          #eval.metric = c('TSS'),
                                          #eval.metric.quality.threshold = quantileThresh,
                                          prob.mean = TRUE,
                                          prob.mean.weight.decay = 'proportional'))


# Get evaluation scores for the Ensemble Modelling stage
emEvalDF <- as.data.frame(get_evaluations(myBiomodEM))
write.csv(emEvalDF, file = paste(getwd(),"/",sp,"/",sp,"_EnsMod_evalDF_AllMetrics.csv",sep=""))


## -------------------------------------------------------------------------------------- ##
## Obtain spatiotemporal projections ----
## -------------------------------------------------------------------------------------- ##

# Models to consider in the ensemble and projection
modelsToUse <- get_kept_models(myBiomodEM, 1)



for(projName in projNames){
  
  cat("## -------------------------------------------------------------------------- ##\n")
  cat(toupper("Making projections for scenario:"), projName,"\n")
  cat("## -------------------------------------------------------------------------- ##\n\n")
  
  # Obtain spatiotemporal projections
  myBiomodProj <- try(BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                        new.env = get(projName),
                                        proj.name = projName, ## Name of the projection from above variable proj.name
                                        selected.models = modelsToUse,
                                        filtered.meth = NULL,
                                        binary.meth = NULL,
                                        compress = 'gzip',
                                        clamping.mask = TRUE,
                                        output.format = '.grd',
                                        do.stack = TRUE))
  
  if(inherits(myBiomodProj, "try-error")){
    sink(file = "log.txt", append=TRUE)
    cat("### ERROR OCCURRED IN SPECIES:",sp,"in BIOMOD_Projection /",projname,"###\n\n")
    #unlink(paste("/home/cibio/myfiles/MountainBryophytesSDM/OUT/MODS/",sp,sep=""), recursive = TRUE)
    sink()
    next
  }
  
  
  # Perform the ensembling of projections
  myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                           binary.meth = c('TSS','ROC','KAPPA'),
                                           EM.output = myBiomodEM,
                                           output.format = '.grd')
  
  # Convert all output raster files to GeoTIFF
  inFolder <- paste(getwd(),"/",sp,"/proj_",projName,sep="")
  outFolder <- paste(inFolder,"/","GeoTIFF", sep="")
  dir.create(outFolder)
  
  if(inherits(myBiomodEF, "try-error")){
    sink(file = "log.txt", append=TRUE)
    cat("### ERROR OCCURRED IN SPECIES:",sp,"in BIOMOD_EnsembleForecasting /",projName,"###\n\n")
    #unlink(paste("/home/cibio/myfiles/MountainBryophytesSDM/OUT/MODS/",sp,sep=""), recursive = TRUE)
    sink()
    next
  }
  
  convertToGeoTIFF(inFolder, outFolder)
  
} 

#save.image(file=paste(sp,"ModObjects-v1.RData",sep="_"))
