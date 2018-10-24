

library(hypervolume)
library(raster)
library(dplyr)
library(magrittr)
library(rgdal)
library(ggplot2)

load("./OUT/HyperVolumeBySpecies-v2.RData")


spDataVars %>% 
  select(-1) %>% 
  cor(method="spearman") %>% 
  round(2) %>% 
  write.csv(file="./RESULTS_2_SHARE/corSpearman_inputDF.csv",row.names=FALSE)

spDataVars %>% 
  select(-1) %>% 
  cor(method="pearson") %>% 
  round(2) %>% 
  write.csv(file="./RESULTS_2_SHARE/corPearson_inputDF.csv",row.names=FALSE)