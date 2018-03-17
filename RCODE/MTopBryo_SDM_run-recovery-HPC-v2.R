

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
load("ESM_Modeling..models.1520298779.out") # GRIALP
#load("ESM_Modeling..models.1520305529.out") # GRICAE
ESM_ModObject <- output

# BIOMOD2 ensemble modelling object (for a specific threshold measure/value)
load("ESM_EnsembleModeling...TSS.0.85.1520298779.out") # GRIALP
#load("ESM_EnsembleModeling...TSS.0.85.1520305529.out") # GRICAE
ESM_EnsembleMod <- output




# Run model projection for a certain raster stack environmental feature space
#
#
#
ESM_Proj <- ecospat.ESM.Projection(ESM.modeling.output = ESM_ModObject,
                                   new.env = a2000,
                                   parallel = TRUE)


ESM_Proj1 <- ecospat.ESM.Projection(ESM.modeling.output = ESM_ModObject,
                                   new.env = a2000,
                                   parallel = FALSE)



# 
# fl_all <- ESM_Proj$pred.biva
# fl <- ESM_Proj$pred.biva[grep(".grd$",ESM_Proj$pred.biva)]
# nblay <- vector(mode = "integer", length = length(fl))
# 
# i<-0
# for(f in fl){
#   i<-i+1
#   rst <- stack(f)
#   nblay[i] <- nlayers(rst)
# }
# 
# 
# 
# r <- raster::stack(fl)
# nl <- nlayers(r)
# 
# failmod <- unlist(ESM_EnsembleMod$failed)
# failmod <- failmod[-grep("none",failmod)]
# 



w <- ESM_EnsembleMod$weights

## proj_a2000.ESM.BIOMOD.1.1520292399
rstFiles <- c()

failedBivCombns <- c()

for(i in ESM_Proj$which.biva){
  
  dirPath <- paste(ESM_ModObject$wd,"/ESM.BIOMOD.",i,"/proj_",ESM_Proj$name.env,".ESM.BIOMOD.",i,".",ESM_Proj$modeling.id,sep="")
  tmpFiles <- list.files(dirPath,pattern=".grd$|.gri$",full.names = TRUE)
  
  rstFiles <- c(rstFiles, tmpFiles)
  
  if(length(tmpFiles)==0){
    failedBivCombns <- c(failedBivCombns,i)
  }
  
  print(i)
  print(tmpFiles)
}



if(length(failedBivCombns) != 0){
 
  new_w <- w[!(names(w) %in% paste(ESM_Proj$models,".ESM.BIOMOD.",failedBivCombns,sep=""))]
  
  ESM_EnsembleMod$weights <- new_w
  
  ESM_Proj$pred.biva <- rstFiles
}




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

