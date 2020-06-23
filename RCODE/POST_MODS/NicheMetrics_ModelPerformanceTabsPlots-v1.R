

library(dplyr)
library(tidyr)
library(ggplot2)

sps <- list.dirs("./OUT/MODS_v2", recursive = FALSE, full.names = FALSE)

for(i in 1:length(sps)){
  
  sp <- sps[i]
  evalDF_tmp <- read.csv(paste("./OUT/MODS_v2/",sp,"/",sp,"_EnsMod_evalDF_AllMetrics.csv",sep=""), 
                         stringsAsFactors = FALSE)[,1:5]
  colnames(evalDF_tmp) <- c("PerfMetric","PerfValue","Cutoff","Sensitivity","Specificity")
  
  if(i==1){
    evalDF <- data.frame(sp=sp,evalDF_tmp)
  }else{
    evalDF <- rbind(evalDF, data.frame(sp=sp,evalDF_tmp))
  }
}

write.csv(evalDF, "./RESULTS_2_SHARE/NicheMetricsOut/ModelPerformance_AllSpecies-v1.csv",row.names = FALSE)

g <- ggplot(evalDF %>% filter(PerfMetric!="KAPPA"), aes(x=PerfMetric, y=PerfValue)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.05, height=0, alpha=0.2, size=2.5) + 
  labs(title = "Model performance",subtitle = "Area under the receiver operating curce (ROC [0,1])\nTrue-skill statistic (TSS [-1,1])") + 
  xlab("Performance metric") + ylab("Performance value")
  

plot(g)

ggsave(plot = g, filename="./RESULTS_2_SHARE/NicheMetricsOut/ModPerformance-v1.png",width=5,height=6)

g1 <- ggplot(evalDF %>% filter(PerfMetric=="TSS")) + 
  geom_boxplot(aes(x="Specificity", y=Specificity)) + 
  geom_boxplot(aes(x="Sensitivity", y=Sensitivity)) + 
  geom_jitter(aes(x="Specificity", y=Specificity), width = 0.05, height=0, alpha=0.2, size=2.5) + 
  geom_jitter(aes(x="Sensitivity", y=Sensitivity), width = 0.05, height=0, alpha=0.2, size=2.5) + 
  #geom_jitter(width = 0.05, height=0, alpha=0.2, size=2.5) + 
  labs(title = "Model performance",subtitle = "Sensitivity (true positive rate) vs. Specificity (true negative rate) \nby TSS cutoff") + 
  xlab("Performance metric") + ylab("Performance value")
  #facet_wrap(.~PerfMetric)


plot(g1)
ggsave(plot = g1, filename="./RESULTS_2_SHARE/NicheMetricsOut/sensitivityVsSpecificity-v1.png", width=5, height=6)


## ----------------------------------------------------------------- ##

library(raster)
library(sf)



pb <- txtProgressBar(1,length(sps),style=3)

for(i in 1:length(sps)){
  
  sp <- sps[i]
  
  ## Current 
  ##
  fpath1 <- paste("./OUT/MODS_v2/",sp,"/proj_current/GeoTIFF/proj_current_",sp,"_ensemble_TSSbin.tif",sep="")
  fpath2 <- paste("./OUT/MODS_v2/",sp,"/proj_current/GeoTIFF/proj_current_",sp,"_ensemble.tif",sep="")
  rstBin <- raster::stack(fpath1)[[1]]
  rstEns <- raster::stack(fpath2)[[1]]/1000
  
  ## Projections of climate change scenarios
  ##
  fpath3 <- paste("./OUT/MODS_v2/",sp,"/proj_he45bi50/GeoTIFF/proj_he45bi50_",sp,"_ensemble_TSSbin.tif",sep="")
  fpath4 <- paste("./OUT/MODS_v2/",sp,"/proj_he85bi50/GeoTIFF/proj_he85bi50_",sp,"_ensemble_TSSbin.tif",sep="")
  fpath5 <- paste("./OUT/MODS_v2/",sp,"/proj_mp45bi50/GeoTIFF/proj_mp45bi50_",sp,"_ensemble_TSSbin.tif",sep="")
  fpath6 <- paste("./OUT/MODS_v2/",sp,"/proj_mp85bi50/GeoTIFF/proj_mp85bi50_",sp,"_ensemble_TSSbin.tif",sep="")
  
  fpath7 <- paste("./OUT/MODS_v2/",sp,"/proj_he45bi50/GeoTIFF/proj_he45bi50_",sp,"_ensemble.tif",sep="")
  fpath8 <- paste("./OUT/MODS_v2/",sp,"/proj_he85bi50/GeoTIFF/proj_he85bi50_",sp,"_ensemble.tif",sep="")
  fpath9 <- paste("./OUT/MODS_v2/",sp,"/proj_mp45bi50/GeoTIFF/proj_mp45bi50_",sp,"_ensemble.tif",sep="")
  fpath10 <- paste("./OUT/MODS_v2/",sp,"/proj_mp85bi50/GeoTIFF/proj_mp85bi50_",sp,"_ensemble.tif",sep="")
  
  rstBin_he45bi50 <- raster::stack(fpath3)[[1]]
  rstBin_he85bi50 <- raster::stack(fpath4)[[1]]
  rstBin_mp45bi50 <- raster::stack(fpath5)[[1]]
  rstBin_mp85bi50 <- raster::stack(fpath6)[[1]]
  
  rstEns_he45bi50 <- raster::stack(fpath7)[[1]]
  rstEns_he85bi50 <- raster::stack(fpath8)[[1]]
  rstEns_mp45bi50 <- raster::stack(fpath9)[[1]]
  rstEns_mp85bi50 <- raster::stack(fpath10)[[1]]
  
  png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/",sp,"_DitribMap-v1.png"),
      res = 300, width=3700, height=1400)
  par(mfcol = c(1,2), mar = c(2,2,2,5))
  plot(rstEns, main = paste(sp,"| Habitat suitability"))
  plot(rstBin)
  dev.off()
  
  png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/",sp,"_DitribMap_he45bi50-v1.png"),
      res = 300, width=3700, height=1400)
  par(mfcol = c(1,2), mar = c(2,2,2,5))
  plot(rstEns_he45bi50, main = paste(sp,"| Habitat suitability (HE/RCP-4.5/2050)"))
  plot(rstBin_he45bi50)
  dev.off()
  
  png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/",sp,"_DitribMap_he85bi50-v1.png"),
      res = 300, width=3700, height=1400)
  par(mfcol = c(1,2), mar = c(2,2,2,5))
  plot(rstEns_he85bi50, main = paste(sp,"| Habitat suitability (HE/RCP-8.5/2050)"))
  plot(rstBin_he85bi50)
  dev.off()
  
  png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/",sp,"_DitribMap_mp45bi50-v1.png"),
      res = 300, width=3700, height=1400)
  par(mfcol = c(1,2), mar = c(2,2,2,5))
  plot(rstEns_mp45bi50, main = paste(sp,"| Habitat suitability (MP/RCP-4.5/2050)"))
  plot(rstBin_mp45bi50)
  dev.off()
  
  png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/",sp,"_DitribMap_mp85bi50-v1.png"),
      res = 300, width=3700, height=1400)
  par(mfcol = c(1,2), mar = c(2,2,2,5))
  plot(rstEns_mp85bi50, main = paste(sp,"| Habitat suitability (MP/RCP-4.5/2050)"))
  plot(rstBin_mp85bi50)
  dev.off()
  
  if(i==1){
    evalDF <- data.frame(sp=sp,evalDF_tmp)
  }else{
    evalDF <- rbind(evalDF, data.frame(sp=sp,evalDF_tmp))
  }
  
  if(i==1){
    
    rstBinSpecies <- rstBin
    
    rstBinSpecies_he45bi50 <- rstBin_he45bi50
    rstBinSpecies_he85bi50 <- rstBin_he85bi50
    rstBinSpecies_mp45bi50 <- rstBin_mp45bi50
    rstBinSpecies_mp85bi50 <- rstBin_mp85bi50
    
  }else{
    
    if(!compareRaster(rstBin,rstBinSpecies,stopiffalse=FALSE, showwarning=FALSE)){
      
      cat("\n\nResampling raster after verifying differences...\n\n")
      
      rstBin          <- resample(rstBin, rstBinSpecies[[1]], method="ngb")
      rstBin_he45bi50 <- resample(rstBin_he45bi50, rstBinSpecies_he45bi50[[1]], method="ngb")
      rstBin_he85bi50 <- resample(rstBin_he85bi50, rstBinSpecies_he85bi50[[1]], method="ngb")
      rstBin_mp45bi50 <- resample(rstBin_mp45bi50, rstBinSpecies_mp45bi50[[1]], method="ngb")
      rstBin_mp85bi50 <- resample(rstBin_mp85bi50, rstBinSpecies_mp85bi50[[1]], method="ngb")
    }
    
    rstBinSpecies <- stack(rstBinSpecies, rstBin)
    rstBinSpecies_he45bi50 <- stack(rstBinSpecies_he45bi50, rstBin_he45bi50)
    rstBinSpecies_he85bi50 <- stack(rstBinSpecies_he85bi50, rstBin_he85bi50)
    rstBinSpecies_mp45bi50 <- stack(rstBinSpecies_mp45bi50, rstBin_mp45bi50)
    rstBinSpecies_mp85bi50 <- stack(rstBinSpecies_mp85bi50, rstBin_mp85bi50)
  }
  setTxtProgressBar(pb,i)
}

## --------------------------------------------------------------------------------------------------------- ##
## Calculate species richness maps

spRichrst <- calc(rstBinSpecies, fun = sum)
writeRaster(spRichrst,filename="./RESULTS_2_SHARE/NicheMetricsOut/RASTER_DATA/speciesRichness_v1.tif")

spRichrst_he45bi50 <- calc(rstBinSpecies_he45bi50, fun = sum)
writeRaster(spRichrst_he45bi50,filename="./RESULTS_2_SHARE/NicheMetricsOut/RASTER_DATA/speciesRichness_he45bi50_v1.tif")

spRichrst_he85bi50 <- calc(rstBinSpecies_he85bi50, fun = sum)
writeRaster(spRichrst_he85bi50,filename="./RESULTS_2_SHARE/NicheMetricsOut/RASTER_DATA/speciesRichness_he85bi50_v1.tif")

spRichrst_mp45bi50 <- calc(rstBinSpecies_mp45bi50, fun = sum)
writeRaster(spRichrst_mp45bi50,filename="./RESULTS_2_SHARE/NicheMetricsOut/RASTER_DATA/speciesRichness_mp45bi50_v1.tif")

spRichrst_mp85bi50 <- calc(rstBinSpecies_mp85bi50, fun = sum)
writeRaster(spRichrst_mp85bi50,filename="./RESULTS_2_SHARE/NicheMetricsOut/RASTER_DATA/speciesRichness_mp85bi50_v1.tif")

## --------------------------------------------------------------------------------------------------------- ##
## Make figures with species richness maps for current and future projections

png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/_spRichnessAllSpecies-v1.png"), 
    res = 300, width=1700, height=1400)
plot(spRichrst,main="Species richness (current)")
dev.off()

png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/_spRichnessAllSpecies_he45bi50-v1.png"), 
    res = 300, width=1700, height=1400)
plot(spRichrst_he45bi50,main="Species richness (HE/RCP-4.5/2050)")
dev.off()

png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/_spRichnessAllSpecies_he85bi50-v1.png"), 
    res = 300, width=1700, height=1400)
plot(spRichrst_he85bi50,main="Species richness (HE/RCP-8.5/2050)")
dev.off()

png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/_spRichnessAllSpecies_mp45bi50-v1.png"), 
    res = 300, width=1700, height=1400)
plot(spRichrst_mp45bi50,main="Species richness (MP/RCP-4.5/2050)")
dev.off()

png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/_spRichnessAllSpecies_mp85bi50-v1.png"), 
    res = 300, width=1700, height=1400)
plot(spRichrst_mp85bi50,main="Species richness (MP/RCP-8.5/2050)")
dev.off()

## --------------------------------------------------------------------------------------------------------- ##
## Calculate species richness difference maps (turnover)

spRichDiff_he45bi50 <- spRichrst_he45bi50 - spRichrst
spRichDiff_he85bi50 <- spRichrst_he85bi50 - spRichrst
spRichDiff_mp45bi50 <- spRichrst_mp45bi50 - spRichrst
spRichDiff_mp85bi50 <- spRichrst_mp85bi50 - spRichrst

writeRaster(spRichDiff_he45bi50,filename="./RESULTS_2_SHARE/NicheMetricsOut/RASTER_DATA/spRichDiff_he45bi50_v1.tif")
writeRaster(spRichDiff_he85bi50,filename="./RESULTS_2_SHARE/NicheMetricsOut/RASTER_DATA/spRichDiff_he85bi50_v1.tif")
writeRaster(spRichDiff_mp45bi50,filename="./RESULTS_2_SHARE/NicheMetricsOut/RASTER_DATA/spRichDiff_mp45bi50_v1.tif")
writeRaster(spRichDiff_mp85bi50,filename="./RESULTS_2_SHARE/NicheMetricsOut/RASTER_DATA/spRichDiff_mp85bi50_v1.tif")


png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/_spRichDiff_he45bi50-v1.png"), 
    res = 300, width=1700, height=1400)
plot(spRichDiff_he45bi50,main="Species richness turnover (HE/RCP-4.5/2050)")
dev.off()

png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/_spRichDiff_he85bi50-v1.png"), 
    res = 300, width=1700, height=1400)
plot(spRichDiff_he85bi50,main="Species richness turnover (HE/RCP-8.5/2050)")
dev.off()

png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/_spRichDiff_mp45bi50-v1.png"), 
    res = 300, width=1700, height=1400)
plot(spRichDiff_mp45bi50,main="Species richness turnover (MP/RCP-4.5/2050)")
dev.off()

png(filename = paste("./RESULTS_2_SHARE/NicheMetricsOut/DISTRIB_MAPS/_spRichDiff_mp85bi50-v1.png"), 
    res = 300, width=1700, height=1400)
plot(spRichDiff_mp85bi50,main="Species richness turnover (MP/RCP-8.5/2050)")
dev.off()

