

library(raster)
library(ecospat)
library(biomod2)
library(dplyr)
library(magrittr)


setwd("D:/MyDocs/temp/colabHH/ModObjects")

spNames <- list.files()

spPath <- list.files(full.names = TRUE)


for(spName in spNames){
  
  varWeights <- list()
  
  
  load(list.files(paste("./",spName, sep=""), pattern="^ESM_Modeling..models.", full.names = TRUE))
  ESM_ModObject <- output
  
  
  load(list.files(paste("./",spName, sep=""), pattern="^ESM_EnsembleModeling...", full.names = TRUE))
  ESM_EnsProjection <- output
  
  w <- ESM_EnsProjection$weights
  
  bivCombValues <- w %>% 
                names %>% 
                strsplit(split = "\\.") %>% 
                lapply(FUN = function(x) x[4]) %>% 
                unlist %>% 
                as.integer
  
  
  spBivCombns <- lapply(ESM_ModObject$mymodels, FUN = function(x,...) x@expl.var.names)
  
  vars <- spBivCombns %>% unlist %>% unique
  
  for(v in vars)
    varWeights[[v]] <- c()
  
  for(i in 1:length(w)){
    
    combVal <- bivCombValues[i]
    varNames <- spBivCombns[[combVal]]
    
    varWeights[[varNames[1]]] <- c(varWeights[[varNames[1]]], w[i])
    varWeights[[varNames[2]]] <- c(varWeights[[varNames[2]]], w[i])   
  }
  
  varWeightsDF <- varWeights %>% as.data.frame %>% stack
  
  boxplot(varWeights)
  
  g <- ggplot(varWeightsDF,aes(x=ind,y=values)) + 
    geom_boxplot(fill="grey", alpha=0.7) + 
    xlab("Variables") + 
    ylab("Variable importance (based on TSS/ESM model weights)") + 
    ggtitle(spName) + 
    theme_light()
  
  plot(g)
  
  ggsave(filename=paste("../VarImp/",spName,"_varImpBoxplot-v1.png",sep=""), plot = g)
  
}



