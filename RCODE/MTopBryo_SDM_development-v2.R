<<<<<<< HEAD


library(raster)
library(ecospat)
library(biomod2)
library(sp)
library(rgdal)
library(dplyr)
library(ggplot2)
library(parallel)
library(doParallel)

## Output directory
outDir <- "D:/MyDocs/Dropbox/Modelling_mountaintop_bryophytes_climate_change/MountainBryophytesSDM/OUT"

#setwd("./OUT")
setwd(outDir)


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



# Loop through all species ---------------------------------------------------------------------------
#
#
#
for(i in 1:length(spCodes)){
  
  # Selected species
  selSpecies <- spCodes[i]
  
  # Filter data to the selected species
  spDataSelSpecies <- as(spData[spData@data$Cod_esp==selSpecies,],"SpatialPoints")
  
  # Perform data formatting
  BIOMOD_Data <- BIOMOD_FormatingData( resp.var = spDataSelSpecies,
                                       expl.var = a2000,
                                       resp.name = selSpecies,
                                       PA.nb.rep = 1, # cannot be >1...
                                       PA.nb.absences = 10000,
                                       PA.strategy = 'random')
  
  # Model options for biomod2
  BIOMOD_ModOptions <- BIOMOD_ModelingOptions(GAM = list(k=2),
                                              RF = list(ntree=250))
  
  # Calibration of simple bivariate models
  ESM_ModObject <- ecospat.ESM.Modeling(data = BIOMOD_Data,
                                        #models = c('GLM','RF','GAM','GBM','FDA','MARS','MAXENT.Tsuruoka'),
                                        models = c('GLM','RF','GAM'),
                                        models.options = BIOMOD_ModOptions,
                                        NbRunEval = 1,
                                        DataSplit = 80,
                                        weighting.score = c("AUC"),
                                        parallel = TRUE) 
  
  
  ESM_ModelsEvaluations <- lapply(ESM_ModObject$mymodels, get_evaluations)
  
  
  # Evaluation and average of simple bivariate models to ESMs
  ESM_EnsembleMod <- ecospat.ESM.EnsembleModeling(ESM_ModObject,
                                                  weighting.score = c("AUC"),
                                                  threshold = 0.9)
  
  write.csv(ESM_EnsembleMod$ESM.evaluations,
            file = paste(selSpecies,"_ESM_EvaluationMetrics_ensProj.csv",sep=""))
  
  
  
  # Loop through all the raster stacks containing different projections
  #
  #
  #
  for(projName in projList){
    
    
    # Projection of simple bivariate models into new feature space -------------
    # Uses eval / parse due to new.env!!!
    #
    
    cat("\n\nPerforming ESM projection for:",projName,"...")
    
    eval(parse(text=
    paste("ESM_Proj <- ecospat.ESM.Projection(ESM.modeling.output = ESM_ModObject,
                                     new.env = ",projName,",
                                     parallel = TRUE)",sep="")
    ))
    
    cat("done.\n")
    
    
    # Projection of calibrated ESMs into new space -------------
    #
    cat("\n\nPerforming ensemble projection for:",projName,"...")
    
    ESM_ensProjection <- ecospat.ESM.EnsembleProjection(ESM.prediction.output = ESM_Proj,
                                                        ESM.EnsembleModeling.output = ESM_EnsembleMod)
    
    cat("done.\n")
    
    
    # Save R objects -------------
    #
    save(BIOMOD_Data, ESM_EnsembleMod, ESM_EnsembleMod, ESM_Proj, ESM_ensProjection, 
         file = paste(selSpecies,"_ESM_Objects_",projName,".RData",sep=""))
    
    
    # Write raster with ensembles -------------
    #
    writeRaster(ESM_ensProjection, 
                filename = paste(selSpecies,"_ESM_ensProj_",projName,".tif",sep=""),
                overwrite=TRUE)
  }

  cat("\n\n||---------- Finished ensemble projection for:",selSpecies,"----------||\n\n")
  
  
  # Go up one level for OUT folder
  #setwd("../")
  if(basename(getwd()) != "OUT")
    setwd(outDir)
  
}
=======


library(raster)
library(ecospat)
library(biomod2)
library(sp)
library(rgdal)
library(dplyr)
library(ggplot2)
library(parallel)
library(doParallel)

## Output directory
outDir <- "D:/MyDocs/Dropbox/Modelling_mountaintop_bryophytes_climate_change/MountainBryophytesSDM/OUT"

#setwd("./OUT")
setwd(outDir)


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
  fl.climVars <- list.files(fld, pattern = ".tif$", full.names = TRUE)
  
  print(fl.climVars)
  
  assign(x = projList[i], value = stack(c(fl.topoVars, fl.climVars)))
  cat("done.\n\n")
}



# Loop through all species ---------------------------------------------------------------------------
#
#
#
for(i in 1:length(spCodes)){
  
  # Selected species
  selSpecies <- spCodes[i]
  
  # Filter data to the selected species
  spDataSelSpecies <- as(spData[spData@data$Cod_esp==selSpecies,],"SpatialPoints")
  
  # Perform data formatting
  BIOMOD_Data <- BIOMOD_FormatingData( resp.var = spDataSelSpecies,
                                       expl.var = a2000,
                                       resp.name = selSpecies,
                                       PA.nb.rep = 1, # cannot be >1...
                                       PA.nb.absences = 10000,
                                       PA.strategy = 'random')
  
  # Model options for biomod2
  BIOMOD_ModOptions <- BIOMOD_ModelingOptions(GAM = list(k=2),
                                              RF = list(ntree=250))
  
  # Calibration of simple bivariate models
  ESM_ModObject <- ecospat.ESM.Modeling(data = BIOMOD_Data,
                                        #models = c('GLM','RF','GAM','GBM','FDA','MARS','MAXENT.Tsuruoka'),
                                        models = c('GLM','RF','GAM'),
                                        models.options = BIOMOD_ModOptions,
                                        NbRunEval = 1,
                                        DataSplit = 80,
                                        weighting.score = c("AUC"),
                                        parallel = TRUE) 
  
  
  ESM_ModelsEvaluations <- lapply(ESM_ModObject$mymodels, get_evaluations)
  
  
  # Evaluation and average of simple bivariate models to ESMs
  ESM_EnsembleMod <- ecospat.ESM.EnsembleModeling(ESM_ModObject,
                                                  weighting.score = c("AUC"),
                                                  threshold = 0.9)
  
  write.csv(ESM_EnsembleMod$ESM.evaluations,
            file = paste(selSpecies,"_ESM_EvaluationMetrics_ensProj.csv",sep=""))
  
  
  
  # Loop through all the raster stacks containing different projections
  #
  #
  #
  for(projName in projList){
    
    
    # Projection of simple bivariate models into new feature space -------------
    # Uses eval / parse due to new.env!!!
    #
    
    cat("\n\nPerforming ESM projection for:",projName,"...")
    
    eval(parse(text=
    paste("ESM_Proj <- ecospat.ESM.Projection(ESM.modeling.output = ESM_ModObject,
                                     new.env = ",projName,",
                                     parallel = TRUE)",sep="")
    ))
    
    cat("done.\n")
    
    
    # Projection of calibrated ESMs into new space -------------
    #
    cat("\n\nPerforming ensemble projection for:",projName,"...")
    
    ESM_ensProjection <- ecospat.ESM.EnsembleProjection(ESM.prediction.output = ESM_Proj,
                                                        ESM.EnsembleModeling.output = ESM_EnsembleMod)
    
    cat("done.\n")
    
    
    # Save R objects -------------
    #
    save(BIOMOD_Data, ESM_EnsembleMod, ESM_EnsembleMod, ESM_Proj, ESM_ensProjection, 
         file = paste(selSpecies,"_ESM_Objects_",projName,".RData",sep=""))
    
    
    # Write raster with ensembles -------------
    #
    writeRaster(ESM_ensProjection, 
                filename = paste(selSpecies,"_ESM_ensProj_",projName,".tif",sep=""),
                overwrite=TRUE)
  }

  cat("\n\n||---------- Finished ensemble projection for:",selSpecies,"----------||\n\n")
  
  
  # Go up one level for OUT folder
  #setwd("../")
  if(basename(getwd()) != "OUT")
    setwd(outDir)
  
}
>>>>>>> 5b50d849ebd42783f0e0894c5fa6594a332de661
