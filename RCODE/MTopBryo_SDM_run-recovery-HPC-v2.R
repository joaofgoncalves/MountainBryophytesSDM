

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
outDir <- "~/myfiles/MountainBryophytesSDM/OUT/ESM.BIOMOD.output_GRIALP"

#setwd("./OUT")
setwd(outDir)

# BIOMOD2 modelling object
load("ESM_Modeling..models.1520298779.out")
ESM_ModObject <- output

# BIOMOD2 ensemble modelling object (for a specific threshold measure/value)
load("ESM_EnsembleModeling...TSS.0.85.1520298779.out")
ESM_EnsembleMod <- output




# Run model projection for a certain raster stack environmental feature space
#
#
#
ESM_Proj <- ecospat.ESM.Projection(ESM.modeling.output = ESM_ModObject,
                                   new.env = a2000,
                                   parallel = FALSE)



ESM_ensProjection <- ecospat.ESM.EnsembleProjection(ESM.prediction.output = ESM_Proj,
                                                    ESM.EnsembleModeling.output = ESM_EnsembleMod)


# Save R objects -------------
#
save(BIOMOD_Data, ESM_EnsembleMod, ESM_EnsembleMod, ESM_Proj, ESM_ensProjection, 
     file = paste(selSpecies,"_ESM_Objects_",projName,".RData",sep=""))


# Write raster with ensembles -------------
#
writeRaster(ESM_ensProjection, 
            filename = paste(selSpecies,"_ESM_ensProj_",projName,".tif",sep=""),
            overwrite=TRUE)

