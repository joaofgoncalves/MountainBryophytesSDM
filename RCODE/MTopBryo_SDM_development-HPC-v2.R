

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

#cl <- makeCluster(25)
#registerDoParallel(cl)


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


# Model options for biomod2
BIOMOD_ModOptions <- BIOMOD_ModelingOptions(GAM = list(k=4),
                                            RF  = list(ntree=250),
                                            GBM = list(n.trees = 1000,
                                                       interaction.depth = 5),
                                            ANN = list(maxit=150)
)



# Recover previous model runs -----------------------------------------------------------
# 
# selSpecies <- "ANDHEI"
# 
# ## Output directory
# outDir <- "~/myfiles/MountainBryophytesSDM/OUT/ESM.BIOMOD.output_ANDHEI"
# 
# #setwd("./OUT")
# setwd(outDir)
# 
# # BIOMOD2 modelling object
# load("ESM_Modeling..models.1520351644.out") # ANDHEI
# ESM_ModObject <- output
# 
# # BIOMOD2 ensemble modelling object (for a specific threshold measure/value)
# load("ESM_EnsembleModeling...TSS.0.85.1520351644.out") # ANDHEI
# ESM_EnsembleMod <- output
# 
# 



# Loop through all species ---------------------------------------------------------------------------
#
#
#
#
#for(i in 1:length(spCodes)){
for(i in 3){
    
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
  
  
  # Calibration of simple bivariate models
  ESM_ModObject <- try(ecospat.ESM.Modeling(data = BIOMOD_Data,
                                        models = c('GLM','CTA','GAM','RF','ANN','MARS'),
                                        models.options = BIOMOD_ModOptions,
                                        NbRunEval = 10,
                                        DataSplit = 80,
                                        weighting.score = "TSS",
                                        parallel = FALSE))
  
  
  #ESM_ModelsEvaluations <- lapply(ESM_ModObject$mymodels, get_evaluations)
  
  if(inherits(ESM_ModObject,"try-error")){
    cat("\n\n!!!Error while performing model calibration for:", selSpecies,"!!!\n\n\n")
    next
  }
    
  
  # Evaluation and average of simple bivariate models to ESMs
  ESM_EnsembleMod <- ecospat.ESM.EnsembleModeling(ESM_ModObject,
                                                  weighting.score = c("TSS"),
                                                  threshold = 0.85)
  
  write.csv(ESM_EnsembleMod$ESM.evaluations,
            file = paste(selSpecies,"_ESM_EvaluationMetrics_ensProj.csv",sep=""))
  
  
  # Loop through all the raster stacks containing different projections --------
  #
  #
  #
  for(projName in projList){
    
    
    # Projection of simple bivariate models into new feature space -------------
    # Uses eval / parse due to new.env!!!
    #
    
    cat("\n\nPerforming ESM projection for:",projName,"...\n")
    
    eval(parse(text=
    paste("ESM_Proj <- try(ecospat.ESM.Projection(ESM.modeling.output = ESM_ModObject,
                                     new.env = ",projName,",
                                     parallel = FALSE))",sep="")
    ))
    
    
    if(inherits(ESM_Proj,"try-error")){
      cat("\n\n!!!Error while performing model spatiotemporal projections for:", selSpecies,"and projection:",projName,"!!!\n\n\n")
      next
    }
    
    cat("done.\n")
    
    
    # Checking files -----------------------------------------------------------
    #
    
    cat("Checking files and weights:\n")
    
    w <- ESM_EnsembleMod$weights
    rstFiles <- c()
    failedBivCombns <- c()
    
    for(i in ESM_Proj$which.biva){
      
      dirPath <- paste(ESM_ModObject$wd,"/ESM.BIOMOD.",i,"/proj_",ESM_Proj$name.env,".ESM.BIOMOD.",i,".",ESM_Proj$modeling.id,sep="")
      tmpFiles <- list.files(dirPath,pattern=".grd$|.gri$",full.names = TRUE)
      
      rstFiles <- c(rstFiles, tmpFiles)
      
      if(length(tmpFiles)==0){
        failedBivCombns <- c(failedBivCombns,i)
      }
      #print(i)
      #print(tmpFiles)
    }
    
    # Replace metadata in objects that failed to obtain projections
    #
    if(length(failedBivCombns) != 0){
      
      # Create the modified ESM/BIOMOD objects
      ESM_Proj_mod <- ESM_Proj
      ESM_EnsembleMod_mod <- ESM_EnsembleMod
      
      # Model codes to remove / not run or finished during projection step
      modsToRemove <- paste(rep(ESM_Proj$models, each=length(failedBivCombns)),".ESM.BIOMOD.",failedBivCombns,sep="")
      
      # New weights generated by removing the unfinished projections
      new_w <- w[!(names(w) %in% modsToRemove)]
      
      # Input values into the mod objects
      ESM_EnsembleMod_mod$weights <- new_w
      ESM_Proj_mod$pred.biva <- rstFiles
      
      cat("Found",length(failedBivCombns),"files lacking for combns:",paste(failedBivCombns,collapse="/"),"...\n")
    }
    
    cat("done!\n\n")
    
    
    # Projection of calibrated ESMs into new space -------------
    #
    cat("\n\nPerforming ensemble projection for:",projName,"...\n")
    
    
    if(length(failedBivCombns) != 0){
      ESM_ensProjection <- try(ecospat.ESM.EnsembleProjection(ESM.prediction.output = ESM_Proj_mod,
                                                              ESM.EnsembleModeling.output = ESM_EnsembleMod_mod))
    }else{
      ESM_ensProjection <- try(ecospat.ESM.EnsembleProjection(ESM.prediction.output = ESM_Proj,
                                                              ESM.EnsembleModeling.output = ESM_EnsembleMod))
    }  
    
    
    if(inherits(ESM_ensProjection,"try-error")){
      cat("\n\n!!!Error while performing ensemble model projections for:", selSpecies,"and projection:",projName,"!!!\n\n\n")
      next
    }
    
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



