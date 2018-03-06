

library(raster, lib.loc = "~/R/powerpc64le-redhat-linux-gnu-library/3.3")
library(earth, lib.loc = "~/R/powerpc64le-redhat-linux-gnu-library/3.3")
library(gbm, lib.loc = "~/R/powerpc64le-redhat-linux-gnu-library/3.3")
library(ecospat, lib.loc = "~/R/powerpc64le-redhat-linux-gnu-library/3.3")
library(biomod2)
library(sp)
library(rgdal)
library(dplyr)
library(ggplot2)
library(parallel)
library(doParallel)

## Output directory
#outDir <- "D:/MyDocs/Dropbox/Modelling_mountaintop_bryophytes_climate_change/MountainBryophytesSDM/OUT"
outDir <- "~/myfiles/MountainBryophytesSDM/OUT"

#setwd("./OUT")
setwd(outDir)


## Start cluster for parallel processing -------------------------------------------------------------

cl <- makeCluster(25)
registerDoParallel(cl)


## Load data for selected bryophyte species ----------------------------------------------------------


spData <- readOGR("../DATA/VECTOR/Bryophyte_dataset", "And_Gri_Rac_PI_mountaintop", 
                  stringsAsFactors = FALSE)

# Get the species codes
spCodes <- unique(spData@data$Cod_esp)


## Spatial data load ---------------------------------------------------------------------------------


## File list for different climate projections dates and GCMs
##

dirList <- list.dirs("../DATA/RASTER/Worldclim_1km")[-1] #Remove the base dir
projList <- c() # Vector for holding the name for each projection


#fl.topoVars <- list.files("../DATA/RASTER/TopoVars", full.names = TRUE, pattern = ".tif$")[c(1,2,3,5)]
fl.topoVars <- list.files("../DATA/RASTER/TopoVars", full.names = TRUE, pattern = ".tif$")[c(1,5)]


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
#
for(i in 4:length(spCodes)){
#for(i in 10){
  
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
  BIOMOD_ModOptions <- BIOMOD_ModelingOptions(GAM = list(k=4),
                                              RF  = list(ntree=250),
                                              GBM = list(n.trees = 1000,
                                                         interaction.depth = 5),
                                              ANN = list(maxit=150)
                                              )
  
  # Calibration of simple bivariate models
  ESM_ModObject <- ecospat.ESM.Modeling(data = BIOMOD_Data,
                                        #models = c('GLM','RF','GAM','GBM','FDA','MARS','MAXENT.Tsuruoka'),
                                        models = c('GLM','CTA','GAM','RF','ANN','MARS'),
                                        models.options = BIOMOD_ModOptions,
                                        NbRunEval = 10,
                                        DataSplit = 80,
                                        weighting.score = "TSS",
                                        parallel = TRUE) 
  
  
  #ESM_ModelsEvaluations <- lapply(ESM_ModObject$mymodels, get_evaluations)
  
  
  # Evaluation and average of simple bivariate models to ESMs
  ESM_EnsembleMod <- ecospat.ESM.EnsembleModeling(ESM_ModObject,
                                                  weighting.score = c("TSS"),
                                                  threshold = 0.85)
  
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
    
    cat("\n\nPerforming ESM projection for:",projName,"...\n")
    
    eval(parse(text=
    paste("ESM_Proj <- ecospat.ESM.Projection(ESM.modeling.output = ESM_ModObject,
                                     new.env = ",projName,",
                                     parallel = TRUE)",sep="")
    ))
    
    cat("done.\n")
    
    
    # Projection of calibrated ESMs into new space -------------
    #
    cat("\n\nPerforming ensemble projection for:",projName,"...\n")
    
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
