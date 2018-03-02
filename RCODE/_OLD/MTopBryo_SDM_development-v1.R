

library(raster)
library(ecospat)
library(biomod2)
library(sp)
library(rgdal)
library(dplyr)
library(ggplot2)
library(parallel)
library(doParallel)


setwd("./OUT")



## Start cluster for parallel processing -------------------------------------------------------------

cl <- makeCluster(5)
registerDoParallel(cl)


## Load data for selected bryophyte species ----------------------------------------------------------


spData <- readOGR("../DATA/VECTOR/Bryophyte_dataset", "And_Gri_Rac_PI_mountaintop", 
                  stringsAsFactors = FALSE)

spCodes <- unique(spData@data$Cod_esp)


## Spatial data load ---------------------------------------------------------------------------------


## File list for different climate projections dates and GCMs
##

dirList <- list.dirs("../DATA/RASTER/Worldclim_1km")[-1]
projList <- c()


fl.topoVars <- list.files("../DATA/RASTER/TopoVars", full.names = TRUE, pattern = ".tif$")[c(1,2,3,5)]


i <- 0
for(fld in dirList){

  cat("-> Loading data for:",fld,"..... ")  
  i<-i+1
  projList[i] <- basename(fld)
  fl.climVars <- list.files(fld,pattern = ".tif$", full.names = TRUE)
  assign(x = projList[i], value = stack(c(fl.topoVars, fl.climVars)))
  cat("done.\n\n")
}




for(i in 1:length(spCodes)){
  
  # Selected species
  selSpecies <- spCodes[i]
  
  # Filter data to the selected species
  spDataSelSpecies <- as(spData[spData@data$Cod_esp==selSpecies,],"SpatialPoints")
  
  # Perform data formatting
  BIOMOD_Data <- BIOMOD_FormatingData( resp.var = spDataSelSpecies,
                                       expl.var = a2000,
                                       resp.name = selSpecies,
                                       #PA.nb.rep = 1, # cannot be >1...
                                       PA.nb.absences = 10000,
                                       PA.strategy = 'random')
  
  # Model options for biomod2
  BIOMOD_ModOptions <- BIOMOD_ModelingOptions(GAM = list(k=2),
                                              RF = list(ntree=250))
  
  # Calibration of simple bivariate models
  ESM_ModObject <- ecospat.ESM.Modeling(data = BIOMOD_Data,
                                        #models = c('GLM','RF','GAM','GBM','FDA','MARS','MAXENT.Tsuruoka'),
                                        models = c('GLM','RF'),
                                        models.options = BIOMOD_ModOptions,
                                        NbRunEval = 1,
                                        DataSplit = 75,
                                        weighting.score = c("AUC"),
                                        parallel = FALSE) 
  
  
  ESM_ModelsEvaluations <- lapply(ESM_ModObject$mymodels, get_evaluations)
  
  
  # Evaluation and average of simple bivariate models to ESMs
  ESM_EnsembleMod <- ecospat.ESM.EnsembleModeling(ESM_ModObject,
                                                  weighting.score = c("AUC"),
                                                  threshold = 0.9975)
  
  # 
  # write.csv(ESM_EnsembleMod$ESM.evaluations,
  #           file = paste("./ESM.BIOMOD.output_",selSpecies,"/ESM_ensProj_",selSpecies,".csv",sep=""))
  # 
  # 
  write.csv(ESM_EnsembleMod$ESM.evaluations,
            file = paste("ESM_ensProj_",selSpecies,".csv",sep=""))
  
  
  for(projName in projList){
    
    
    # Projection of simple bivariate models into new space 
    
    eval(parse(text=
    paste("ESM_Proj <- ecospat.ESM.Projection(ESM.modeling.output = ESM_ModObject,
                                     new.env = ",projName,",
                                     parallel = FALSE)",sep="")
    ))
    
    # Projection of calibrated ESMs into new space 
    ESM_ensProjection <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=ESM_Proj,
                                                        ESM.EnsembleModeling.output=ESM_EnsembleMod)
    
    # # Save R objects 
    # save(BIOMOD_Data, ESM_object, ESM_EnsembleMod, ESM_Proj, ESM_ensProjection, 
    #      file = paste("./ESM.BIOMOD.output_",selSpecies,"/ESM_Objects_",projName,"_",selSpecies,".RData",sep=""))
    # 
    # # Write raster with ensembles
    # writeRaster(ESM_ensProjection, 
    #             filename = paste("./ESM.BIOMOD.output_",selSpecies,"/ESM_ensProj_",projName,"_",selSpecies,".tif",sep=""))
    # 
    # 
    # 
    
    # Save R objects 
    save(BIOMOD_Data, ESM_EnsembleMod, ESM_EnsembleMod, ESM_Proj, ESM_ensProjection, 
         file = paste("ESM_Objects_",projName,"_",selSpecies,".RData",sep=""))
    
    # Write raster with ensembles
    writeRaster(ESM_ensProjection, 
                filename = paste("ESM_ensProj_",projName,"_",selSpecies,".tif",sep=""),
                overwrite=TRUE)
    
    
    
  }
  
  
  # Go up one level for OUT folder
  #setwd("../")
  setwd("../../")
}




