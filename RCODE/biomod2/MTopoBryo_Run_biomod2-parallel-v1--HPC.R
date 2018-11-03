
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
#library(foreach)
#library(parallel)
#library(doParallel)

# cl <- makeCluster(3)
# registerDoParallel(cl)


#setwd("D:/MyDocs/Dropbox/Modelling_mountaintop_bryophytes_climate_change/MountainBryophytesSDM")
setwd("~/myfiles/MountainBryophytesSDM")

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
climVarsPaths <- list.files("./DATA/RASTER/Clim/a2000", pattern=".tif$", full.names = TRUE)
topoVarsPaths <- list.files("./DATA/RASTER/TopoVars", pattern=".tif$", full.names = TRUE)[c(2,5)]
soilVarsPaths <- list.files("./DATA/RASTER/Soil", pattern=".tif$", full.names = TRUE)[c(8)]

# Load raster stack with climatic and topographic data
current <- stack(c(climVarsPaths,topoVarsPaths,soilVarsPaths))
names(current) <- c(paste("BIO",c(11,17,19),sep=""),"ASPBR","TOPRI","SOIPH")

# Directories 
#climProjDirs <- list.dirs("./DATA/RASTER/WorldClim_1km", recursive = FALSE)[-c(1,seq(3,9,by=2))]
climProjDirs <- list.dirs("./DATA/RASTER/Clim", recursive = FALSE)[-c(1,seq(3,9,by=2))]

climProjNames <- basename(climProjDirs)


for(projDir in climProjDirs){
  
  climProjName <- basename(projDir)
  print(climProjName)
  
  assign(climProjName, stack(list.files(projDir, pattern = ".tif$", full.names = TRUE)))
  
  assign(climProjName, stack(projectRaster(get(climProjName),current, method = "ngb"), 
                             stack(topoVarsPaths,soilVarsPaths)))
  
  assign(climProjName,`names<-`(get(climProjName), 
                                c(paste("BIO",c(11,17,19),sep=""),"ASPBR","TOPRI","SOIPH")))
  
}


## -------------------------------------------------------------------------------------- ##
## PARAMETERS ----
## -------------------------------------------------------------------------------------- ##

projNames <- c("current", climProjNames)

# Model hyperparameters
# GAM: changes k=4 to avoid overly complex models
#
myBiomodOptions <- biomod2::BIOMOD_ModelingOptions(GAM = list(k = 4),
                                                   MAXENT.Phillips = list(threshold=FALSE,
                                                                          hinge=FALSE,
                                                                          #path_to_maxent.jar="D:/MyDocs/temp"),
                                                                          path_to_maxent.jar="~/myfiles/MountainBryophytesSDM"),
                                                   GBM = list(n.trees = 1500,
                                                              n.cores = 1))

## -------------------------------------------------------------------------------------- ##
## Run biomod2 ----
## -------------------------------------------------------------------------------------- ##


for(i in 5:length(allSpNames)){
#foreach(i = 1:length(allSpNames), .verbose  = TRUE, 
  # .packages = c("raster","biomod2","sp","dplyr","magrittr","biomod2plus"), 
  # .export   = c(projNames,"projNames","allSpNames","myBiomodOptions","spBryoDF")) %dopar% {
  # 
  
  setwd("~/myfiles/MountainBryophytesSDM/OUT/MODS")                                  
  
  sp <- as.character(allSpNames)[i]                                  
  #sp <- abbrevNames(spName)
  
  spRecordPoints <- SpatialPoints(spBryoDF %>% 
                                    filter(Cod_esp == sp) %>% 
                                    dplyr::select(lon, lat),
                                  proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  
  #spRecordPoints <- getUniqueRasterXYCoords(current, spRecordPoints, spatial = TRUE)
  
  # Number of training points
  Npresences <- length(spRecordPoints)
  
  # Set up the biomod data object for calibration
  # Number of PA sets = 10
  # Number of PA's per set = Number of presences
  # PA selection strategy = random
  #
  myBiomodData <- BIOMOD_FormatingData(resp.var = spRecordPoints,
                                       expl.var = current,
                                       resp.name = sp,
                                       PA.nb.rep = 3, # Number of pseudo-absences sets
                                       #PA.nb.absences = ifelse(Npresences < 500, Npresences*10, Npresences), # Nr of pseudo-absences
                                       PA.nb.absences = biomod2plus:::powerFunPAnumberCalc(Npresences),
                                       PA.strategy = 'random') # PA generation method
  
  ## -------------------------------------------------------------------------------------- ##
  ## Calibrate models ----
  ## -------------------------------------------------------------------------------------- ##
  
  
  myBiomodModelOut <- try(BIOMOD_Modeling(
    data = myBiomodData, # Input data
    models = c('GLM',
               'GBM',
               'GAM','CTA','ANN',
               'FDA','MARS','RF',
               'MAXENT.Phillips', 
               'MAXENT.Tsuruoka'), # Models to run
    models.options = myBiomodOptions,
    NbRunEval = 10, # Number of Evaluation runs
    DataSplit = biomod2plus:::adjustPercTrain(myBiomodData, 
                                              percTrain = 90, maxTrain = 2000),
    Prevalence = 0.5, # Prevalence between 0 and 1
    VarImport = 3, # Nr of rounds to evaluate variables
    models.eval.meth = c('TSS','ROC','KAPPA'), # Evaluation metrics
    SaveObj = TRUE,
    rescal.all.models = FALSE,
    do.full.models = TRUE) # Model with all data?
  )
  
  
  if(inherits(myBiomodModelOut, "try-error")){
    sink(file = "log.txt", append=TRUE)
    cat("### ERROR OCCURRED IN SPECIES:",sp,"in BIOMOD_Modeling ###\n\n")
    #unlink(paste("~/myfiles/MountainBryophytesSDM/OUT/MODS/",sp,sep=""), recursive = TRUE)
    sink()
    next
  }
  
  
  # Get model evaluation values
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  
  # Print ROC scores
  print(myBiomodModelEval["ROC","Testing.data",,,])
  print(myBiomodModelEval["TSS","Testing.data",,,])
  
  # Get boxplot stats
  print(fivenum(as.numeric(myBiomodModelEval["ROC","Testing.data",,,])))
  print(fivenum(as.numeric(myBiomodModelEval["TSS","Testing.data",,,])))
  
  # Save evaluation metrics from the arrays
  evalDF.ROC <- as.data.frame(myBiomodModelEval["ROC","Testing.data",,,])
  evalDF.TSS <- as.data.frame(myBiomodModelEval["TSS","Testing.data",,,])
  evalDF.KAPPA <- as.data.frame(myBiomodModelEval["KAPPA","Testing.data",,,])
  
  write.csv(evalDF.ROC, file = paste(getwd(),"/",sp,"/",sp,"_evalDF_ROC.csv",sep=""))
  write.csv(evalDF.TSS, file = paste(getwd(),"/",sp,"/",sp,"_evalDF_TSS.csv",sep=""))
  write.csv(evalDF.KAPPA, file = paste(getwd(),"/",sp,"/",sp,"_evalDF_KAPPA.csv",sep=""))
  
  
  # Calculate variable importance across all PA sets, eval rouns and algorithms 
  varImportance <- get_variables_importance(myBiomodModelOut)
  varImportanceByVariableAVG <- apply(varImportance,1,mean, na.rm = TRUE)
  varImportanceByVariableSTD <- apply(varImportance,1,sd, na.rm = TRUE)
  
  vimpDF <- data.frame(cnames=names(varImportanceByVariableAVG),
                       vimpAVG = varImportanceByVariableAVG, 
                       varImpSTD=varImportanceByVariableSTD) %>% 
    arrange(desc(vimpAVG))
  
  print(vimpDF)
  write.csv(vimpDF, file = paste(getwd(),"/",sp,"/",sp,"_varImportance.csv",sep=""))
  
  try(
    varImportancePlot(myBiomodModelOut,"all",outputFolder = paste("./",sp,sep=""),
                      filename = paste(sp,"_VarImpPlot-all.png",sep=""), plot=FALSE, save=TRUE)
  )
  
  ## Make response plots for all variables ----
  
  #allModAlgos <- c("GBM", "GAM", "RF")
  # allModAlgos <- names(sort(apply(myBiomodModelEval["TSS","Testing.data",,,],1,mean),
  #                           decreasing = TRUE))[1:3]
  # 
  # for(modAlgo in allModAlgos){
  #   
  #   cat("Processing algorithm:",modAlgo,".....\n\n")
  #   
  #   # Load required models
  #   myMods <- BIOMOD_LoadModels(myBiomodModelOut, models=modAlgo)
  #   
  #   respPlotDir <- paste("./",sp,sep="")
  #   if(!dir.exists(respPlotDir)){
  #     dir.create(respPlotDir)
  #   }
  #   
  #   # Make response plots for all variables
  #   try(responsePlots(myBiomodModelOut, Data = get_formal_data(myBiomodModelOut,'expl.var'), modelsToUse = myMods, 
  #                     showVars = "all", fixedVarMetric = 'mean', plotStdErr = TRUE,
  #                     addMarginalPlot = FALSE, marginPlotType = "histogram",
  #                     filePrefix = paste(sp,"_RespPlot_",sep=""), fileSuffix = paste("_",modAlgo,sep=""),
  #                     outFolder = respPlotDir, height = 4, 
  #                     width = 4, plot = FALSE, save = TRUE))
  #   
  #   rm(list = myMods)
  #   cat(" done.\n\n")
  # }
  # 
  
  ## -------------------------------------------------------------------------------------- ##
  ## Perform ensemble modelling ----
  ## -------------------------------------------------------------------------------------- ##
  
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
  
  if(inherits(myBiomodEM, "try-error")){
    sink(file = "log.txt", append=TRUE)
    cat("### ERROR OCCURRED IN SPECIES:",sp,"in BIOMOD_EnsembleModeling ###\n\n")
    #unlink(paste("~/myfiles/MountainBryophytesSDM/OUT/MODS/",sp,sep=""), recursive = TRUE)
    sink()
    next
  }
  
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
      #unlink(paste("~/myfiles/MountainBryophytesSDM/OUT/MODS/",sp,sep=""), recursive = TRUE)
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
      #unlink(paste("~/myfiles/MountainBryophytesSDM/OUT/MODS/",sp,sep=""), recursive = TRUE)
      sink()
      next
    }
    
    convertToGeoTIFF(inFolder, outFolder)
    
  } 
  
  save.image(file=paste(sp,"ModObjects-v1.RData",sep="_"))
}






