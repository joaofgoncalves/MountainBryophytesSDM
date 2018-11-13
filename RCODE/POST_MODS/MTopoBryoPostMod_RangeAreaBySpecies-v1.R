
library(raster)
library(biomod2)
library(ggplot2)
library(dplyr)

fl <- list.files("/mnt/dados/datasets/jg/OUT/MODS", pattern=".RData$", full.names = TRUE)

spNames <- list.dirs("/mnt/dados/datasets/jg/OUT/MODS", full.names = FALSE, recursive = FALSE)

dirs <- list.dirs("/mnt/dados/datasets/jg/OUT/MODS", full.names = TRUE, recursive = FALSE)

projNames <- c("current","he45bi50","he85bi50","mp45bi50","mp85bi50")


for(projName in projNames){
  
  cat("\n\n-> Working with projection:",projName,"......\n\n")
  pb <- txtProgressBar(min=1,max=length(spNames),style=3)
  
  for(i in 1:length(spNames)){
    
    rstProj <- paste(dirs[i],"/proj_",projName,"/GeoTIFF/proj_",projName,"_",spNames[i],"_ensemble_TSSbin.tif", sep="")
    
    rstStackProj <- stack(rstProj)
    
    rstDF <- na.omit(values(rstStackProj))
    
    dfTmp <- data.frame(spNames[i], t(apply(rstDF[,1:3],2,sum)), nrow(rstDF))
    colnames(dfTmp) <- c("spNames", "suitArea_TSS", "suitArea_ROC", "suitArea_KAPPA", "totalArea_npix")
    
    if(i==1){
      DF_ProjArea <- dfTmp
    }else{
      DF_ProjArea <- rbind(DF_ProjArea, dfTmp)
    }
    
    setTxtProgressBar(pb,i)
  }
  
  write.csv(DF_ProjArea, file = paste("ProjectedAreaBySpecies_proj_",projName,"-v1.csv",sep=""), row.names = FALSE)
  rm(DF_ProjArea)
  
}

