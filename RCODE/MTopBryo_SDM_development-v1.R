

library(raster)
library(ecospat)
library(biomod2)
library(sp)
library(rgdal)
library(dplyr)
library(ggplot2)

setwd("./OUT")

spData <- readOGR("./DATA/VECTOR/Bryophyte_dataset", "And_Gri_Rac_PI_mountaintop", 
                  stringsAsFactors = FALSE)

spCodes <- unique(spData@data$Cod_esp)



r <- stack(c("../DATA/RASTER/Worldclim_1km/a2000/bio_03.tif",
             "D:/MyDocs/temp/colabHH/Topo1km/Aspect_Beer_1km.tif"))

fl.Current <- list.files("../DATA/RASTER/Worldclim_1km/a2000", pattern=".tif$", full.names = TRUE)

rstStack.Current <- stack(fl.Current)




selSpecies <- spCodes[1]

# Filter data to the selected species
spDataSelSpecies <- as(spData[spData@data$Cod_esp==selSpecies,],"SpatialPoints")


BIOMOD_Data <- BIOMOD_FormatingData( resp.var = spDataSelSpecies,
                                      expl.var = rstStack.Current,
                                      resp.name = selSpecies,
                                      PA.nb.rep = 1,
                                      PA.nb.absences = 10000,
                                      PA.strategy = 'random')

BIOMOD_ModOptions <- BIOMOD_ModelingOptions(GAM=list(k=2),
                                         RF=list(ntree=250))

### Calibration of simple bivariate models
ESM_object <- ecospat.ESM.Modeling( data=BIOMOD_Data,
                                models=c('GLM','RF'),
                                models.options=BIOMOD_ModOptions,
                                NbRunEval=1,
                                DataSplit=80,
                                weighting.score=c("AUC"),
                                parallel=FALSE) 


### Evaluation and average of simple bivariate models to ESMs
ESM_EnsembleMod <- ecospat.ESM.EnsembleModeling(ESM_object,
                                                weighting.score=c("AUC"),
                                                threshold=0.95)

### Projection of simple bivariate models into new space 
ESM_Proj<-ecospat.ESM.Projection(ESM.modeling.output=ESM_object,
                                            new.env=rstStack.Current)

### Projection of calibrated ESMs into new space 
ESM_ensemProjection <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=ESM_Proj,
                                                        ESM.EnsembleModeling.output=ESM_EnsembleMod)


