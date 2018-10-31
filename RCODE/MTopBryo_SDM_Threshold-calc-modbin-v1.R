
library(raster)
library(dplyr)
library(ggplot2)
library(biomod2)
library(ecospat)
library(sp)
library(rgdal)

## ------------------------------------------------------------------------------------------------- ##

outDir <- "D:/MyDocs/Dropbox/Modelling_mountaintop_bryophytes_climate_change/MountainBryophytesSDM/OUT"

setwd(outDir)


## Load data for selected bryophyte species 
spData <- readOGR("../DATA/VECTOR/Bryophyte_dataset", "And_Gri_Rac_PI_mountaintop", 
                  stringsAsFactors = FALSE)

# Get the species codes
spCodes <- unique(spData@data$Cod_esp)

# List of projections
projList <- basename(list.dirs("../DATA/RASTER/Worldclim_1km")[-1]) #Remove the base dir


## ------------------------------------------------------------------------------------------------- ##


maxTSS_Thresh <- function(x){
  
  out <- vector(mode="numeric",length = ncol(x) - 2)
  
  for(i in 3:ncol(x)){
    maxTSS <- suppressWarnings(ecospat.max.tss(DF[,i],DF[,2]))
    out[i-2] <- round(maxTSS[[1]][which.max(maxTSS[[1]][,1]), 2] * 1000)
  }
  return(out)
}


## ------------------------------------------------------------------------------------------------- ##


setwd("D:/MyDocs/temp/colabHH/Mtop_bryophytes")


for(spCode in spCodes){
  
  cat("## ---- ",spCode,"| Loading data and calculating thresholds for calibration data/period ... ")
  
  currentProj_Filename <- paste(spCode,"_ESM_Objects_a2000.RData",sep="")
  
  load(currentProj_Filename)
  
  # Calculate average performance metrics
  ESMevalMetrics_avg <- ESM_EnsembleMod$ESM.evaluations %>%
                            group_by(technique) %>% 
                            summarise_all(mean, na.rm=TRUE) %>% 
                            as.data.frame %>% 
                            .[,-c(2,16,3,8:11)]

  # Calibration points
  pts <- SpatialPoints(BIOMOD_Data@coord, 
                       proj4string = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  # Extract values to calibration points
  ensProjValues <- extract(ESM_ensProjection, pts) / 1000
  
  # Assemble data to calculate thresholds (Observed+PA _vs_ Predicted)
  obs <- BIOMOD_Data@data.species
  DF <- data.frame(ID = 1:length(BIOMOD_Data@data.species),
                   Obs = replace(obs, is.na(obs), 0),
                   ensProjValues)
  
  # If some technique failed remove it from DF (no thresholds calculated...)
  colToRemove <- apply(DF,2,function(x) !all(is.na(x)))
  DF <- DF[, colToRemove]
    
  # Calculate thresholds
  OptThresh <- suppressWarnings(PresenceAbsence::optimal.thresholds(DF, threshold = seq(0, 1, length.out = 1000)))
  OptThresh <- data.frame(method=as.character(OptThresh[,1]),round(OptThresh[,-1] * 1000), stringsAsFactors = FALSE)
  
  # Add maximum TSS threshold
  OptThresh[13, 1] <- "maxTSS"
  OptThresh[13, 2:ncol(OptThresh)] <- maxTSS_Thresh(DF)
  
  # Write threshold data and evaluation metrics
  write.csv(OptThresh, file = paste("./_bin/tables/",spCode,"_OptimalThresholds.csv",sep=""), row.names = FALSE)
  write.csv(ESMevalMetrics_avg, file = paste("./_bin/tables/",spCode,"_ESM_evaluationMetrics_avg.csv",sep=""), row.names = FALSE)
  
  
  cat("done ---- ##\n\n")
  
  ## Iterate by projection stack ----------------------------
  ##
  ##
  
  colToRemove <- colToRemove[-c(1:2)]
  
  for(proj in projList){
    
    cat("-> Working on:",spCode,"| Projection:",proj,".......")
    
    currentProj_Filename <- paste(spCode,"_ESM_Objects_",proj,".RData",sep="")
    
    load(currentProj_Filename)
    
    # Make objects to put binarized predictions
    ESM_ensProjection_bin_05Tresh <- ESM_ensProjection
    ESM_ensProjection_bin_MinROCdist <- ESM_ensProjection
    ESM_ensProjection_bin_maxTSS <- ESM_ensProjection
    
    nl <- nlayers(ESM_ensProjection)
    
    i<-0
    for(j in 1:nl){
      
      if(!colToRemove[j]){
        next
      }else{
        
        i <- i + 1
        
        #cat("Binarizing layer....",i,"....")
        
        MinROCdist_Thresh <- OptThresh[9, i+1] # Min ROC p(0,1) distance
        maxTSS_thresh <- OptThresh[13, i+1] # Max TSS thresh
        
        # Create binary maps
        ESM_ensProjection_bin_05Tresh[[j]] <- ESM_ensProjection_bin_05Tresh[[j]] >= 500
        ESM_ensProjection_bin_MinROCdist[[j]] <- ESM_ensProjection_bin_MinROCdist[[j]] >= MinROCdist_Thresh
        ESM_ensProjection_bin_maxTSS[[j]] <- ESM_ensProjection_bin_maxTSS[[j]] >= maxTSS_thresh
      }
      

      
      #cat("done.\n\n")
    }
    
    # 
    # writeRaster(ESM_ensProjection_bin_05Tresh, 
    #             filename = paste("./_bin/",spCode,"_",proj,"_ESM_ensProj_bin_05Thresh.tif",sep=""))
    # 
    # writeRaster(ESM_ensProjection_bin_MinROCdist, 
    #             filename = paste("./_bin/",spCode,"_",proj,"_ESM_ensProj_bin_minROCdist.tif",sep=""))
    # 
    # writeRaster(ESM_ensProjection_bin_maxTSS, 
    #             filename = paste("./_bin/",spCode,"_",proj,"_ESM_ensProj_bin_maxTSS.tif",sep=""))
    # 
    allBinThreshs <- stack(ESM_ensProjection_bin_05Tresh[[nl]],
                           ESM_ensProjection_bin_MinROCdist[[nl]],
                           ESM_ensProjection_bin_maxTSS[[nl]])
    
    writeRaster(allBinThreshs, 
                filename = paste("./_bin/",spCode,"_",proj,"_ESM_ensProj_bin_EFonly_allThreshs.tif",sep=""))
    
    cat("done.\n\n")
    
  }

}






